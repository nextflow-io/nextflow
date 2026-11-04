/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.pipeline

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.concurrent.ConcurrentHashMap
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.module.spi.RemoteModuleResolver

/**
 * Resolves, installs, lists and removes registry pipelines stored under {@code <baseDir>/pipelines}
 */
@Slf4j
@CompileStatic
class PipelineManager {

    static final Pattern NAME_PATTERN = ~/^[a-z0-9][a-z0-9._-]*\/[a-z][a-z0-9._-]*$/

    private static final String DEFAULT_WORKFLOW = 'main.nf'

    /** An installed pipeline directory with its {@code .pipeline-info} */
    static record InstalledPipeline(Path dir, PipelineInfo info) {}

    /** Outcome of {@link #install}: {@code alreadyInstalled} means nothing was fetched */
    static record InstallResult(PipelineInfo info, boolean alreadyInstalled) {}

    // ponytail: never evicted, fine for a CLI/run lifetime; add a TTL if a long-lived process needs it
    private static final Map<Path,Boolean> MODIFIED = new ConcurrentHashMap<>()

    private final Path baseDir
    private final RegistryClient client
    private final PipelineGitFetcher fetcher

    PipelineManager(Path baseDir, RegistryClient client, PipelineGitFetcher fetcher = new PipelineGitFetcher()) {
        this.baseDir = baseDir
        this.client = client
        this.fetcher = fetcher
    }

    private static boolean isValidName(String name) {
        return name && NAME_PATTERN.matcher(name).matches()
    }

    Path pipelineDir(String name) {
        if( !isValidName(name) )
            throw new AbortOperationException("Invalid pipeline name '$name' -- expected <scope>/<name>, e.g. nf-core/rnaseq")
        return baseDir.resolve(RemoteModuleResolver.PIPELINES_DIR).resolve(name)
    }

    /** Return the main script of a pipeline, installing the latest release when missing */
    Path resolve(String name) {
        final dir = pipelineDir(name)
        final info = PipelineInfo.load(dir)
        if( info == null ) {
            if( Files.exists(dir) ) {
                log.warn1 "Pipeline ${name} has no registry origin (.pipeline-info missing)"
                return dir.resolve(DEFAULT_WORKFLOW)
            }
            return workflowFile(dir, name, install(name, null, false).info().workflowPath)
        }
        if( isModified(info, dir) )
            log.warn1 "Pipeline ${name} has local modifications (checksum mismatch)"
        return workflowFile(dir, name, info.workflowPath ?: DEFAULT_WORKFLOW)
    }

    // the whole-tree checksum is costly and an include resolves the same pipeline several times per run
    private boolean isModified(PipelineInfo info, Path dir) {
        final key = realPath(dir)
        Boolean modified = MODIFIED.get(key)
        if( modified == null ) {
            modified = computeModified(info, dir)
            MODIFIED.put(key, modified)
        }
        return modified
    }

    protected boolean computeModified(PipelineInfo info, Path dir) {
        return info.isModified(dir)
    }

    private static Path realPath(Path dir) {
        return Files.exists(dir) ? dir.toRealPath() : dir.toAbsolutePath().normalize()
    }

    // the workflow path comes from the registry or .pipeline-info, so it must not escape the pipeline dir
    private static Path workflowFile(Path dir, String name, String workflowPath) {
        final file = dir.resolve(workflowPath).normalize()
        if( !file.startsWith(dir) )
            throw new AbortOperationException("Invalid workflow path '${workflowPath}' for pipeline ${name} -- it must stay inside the pipeline directory")
        return file
    }

    InstallResult install(String name, String version, boolean force) {
        final dir = pipelineDir(name)
        final pipeline = client.getPipeline(name)
        final release = version ? client.getPipelineRelease(name, version) : pipeline.latest
        if( release == null )
            throw new AbortOperationException("Pipeline ${name} has no published release")
        if( pipeline.revoked || release.status == PipelineReleaseStatus.REVOKED )
            throw new AbortOperationException("Pipeline ${name}@${release.version} is revoked")

        final workflowPath = pipeline.workflowPath ?: DEFAULT_WORKFLOW
        workflowFile(dir, name, workflowPath)

        final exists = Files.exists(dir)
        final existing = exists ? PipelineInfo.load(dir) : null
        final modified = exists && (existing == null || existing.isModified(dir))
        if( !modified && existing?.version == release.version )
            return new InstallResult(existing, true)
        if( modified && !force )
            throw new AbortOperationException("Pipeline ${name} has local modifications -- use -force to overwrite")

        log.info "Installing pipeline ${name}@${release.version}..."
        final info = new PipelineInfo(
            name: name,
            version: release.version,
            workflowPath: workflowPath,
            source: release.artifact ? 'artifact' : 'git',
            repositoryUrl: pipeline.repositoryUrl,
            revision: release.gitRevision,
            commitId: release.commitId,
            registryUrl: client.primaryUrl(),
            artifactUrl: release.artifact?.url,
            artifactChecksum: release.artifact?.checksum)
        fetchInto(dir, info)
        log.info "Pipeline ${name}@${release.version} installed in ${dir.toAbsolutePath()}"
        return new InstallResult(info, false)
    }

    /** List installed pipelines, skipping in-flight '.tmp-' dirs; entries whose recorded name is invalid are skipped with a warning */
    List<InstalledPipeline> list() {
        final root = baseDir.resolve(RemoteModuleResolver.PIPELINES_DIR)
        final result = new ArrayList<InstalledPipeline>()
        if( !Files.isDirectory(root) )
            return result
        final dirs = Files.walk(root, 2).withCloseable { stream ->
            stream
                .filter { Path p -> root.relativize(p).nameCount == 2 && !p.fileName.toString().startsWith('.') && Files.exists(p.resolve(RemoteModuleResolver.PIPELINE_INFO_FILE)) }
                .sorted()
                .toList()
        } as List<Path>
        for( final dir : dirs ) {
            final info = PipelineInfo.load(dir)
            if( info.name != null && !isValidName(info.name) ) {
                log.warn "Ignoring pipeline in ${dir} -- invalid name '${info.name}' in .pipeline-info"
                continue
            }
            result.add(new InstalledPipeline(dir, info))
        }
        return result
    }

    boolean remove(String name, boolean force) {
        final dir = pipelineDir(name)
        if( !Files.exists(dir) )
            return false
        final info = PipelineInfo.load(dir)
        if( !force && info == null )
            throw new AbortOperationException("Pipeline ${name} has no registry origin (.pipeline-info missing) -- use -force to remove it")
        if( !force && info.isModified(dir) )
            throw new AbortOperationException("Pipeline ${name} has local modifications -- use -force to remove it")
        MODIFIED.remove(realPath(dir))
        FileHelper.deletePath(dir)
        return true
    }

    private void fetchInto(Path dir, PipelineInfo info) {
        Files.createDirectories(dir.parent)
        final tmp = Files.createTempDirectory(dir.parent, ".tmp-${dir.fileName}-")
        try {
            if( info.source == 'artifact' ) {
                // the archive lives in tmp only until extracted, so it never enters the checksum
                final bundle = tmp.resolve('.bundle.tar.gz')
                client.downloadPipelineRelease(info.name, info.version, bundle)
                PipelineBundle.extract(bundle, tmp)
                Files.delete(bundle)
            }
            else {
                fetcher.fetch(info.repositoryUrl, info.revision, info.commitId, tmp)
            }
            info.checksum = PipelineInfo.checksum(tmp)
            info.save(tmp)
            MODIFIED.remove(realPath(dir))
            FileHelper.deletePath(dir)
            Files.move(tmp, dir, StandardCopyOption.ATOMIC_MOVE)
        }
        catch( Throwable e ) {
            try {
                FileHelper.deletePath(tmp)
            }
            catch( IOException ex ) {
                log.warn "Cannot delete temporary directory ${tmp}: ${ex.message}"
            }
            if( e instanceof RegistryException )
                throw new AbortOperationException("Cannot download pipeline ${info.name}@${info.version} -- ${e.message}", e)
            throw e
        }
    }
}
