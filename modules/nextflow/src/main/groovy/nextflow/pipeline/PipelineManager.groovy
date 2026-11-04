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
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import nextflow.exception.AbortOperationException
import nextflow.module.spi.RemoteModuleResolver

/**
 * Resolves, installs, lists and removes registry pipelines stored under {@code <baseDir>/pipelines}
 */
@Slf4j
@CompileStatic
class PipelineManager {

    private static final Pattern NAME_PATTERN = ~/^[a-z0-9][a-z0-9._-]*\/[a-z][a-z0-9._-]*$/

    private static final String DEFAULT_WORKFLOW = 'main.nf'

    private final Path baseDir
    private final RegistryClient client
    private final PipelineGitFetcher fetcher

    PipelineManager(Path baseDir, RegistryClient client, PipelineGitFetcher fetcher = new PipelineGitFetcher()) {
        this.baseDir = baseDir
        this.client = client
        this.fetcher = fetcher
    }

    Path pipelineDir(String name) {
        if( !name || !NAME_PATTERN.matcher(name).matches() )
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
            return workflowFile(dir, name, install(name, null, false).workflowPath)
        }
        if( info.isModified(dir) )
            log.warn1 "Pipeline ${name} has local modifications (checksum mismatch)"
        return workflowFile(dir, name, info.workflowPath ?: DEFAULT_WORKFLOW)
    }

    // the workflow path comes from the registry or .pipeline-info, so it must not escape the pipeline dir
    private static Path workflowFile(Path dir, String name, String workflowPath) {
        final file = dir.resolve(workflowPath).normalize()
        if( !file.startsWith(dir) )
            throw new AbortOperationException("Invalid workflow path '${workflowPath}' for pipeline ${name} -- it must stay inside the pipeline directory")
        return file
    }

    PipelineInfo install(String name, String version, boolean force) {
        final dir = pipelineDir(name)
        final pipeline = client.getPipeline(name)
        final release = version ? client.getPipelineRelease(name, version) : pipeline.latest
        if( release == null )
            throw new AbortOperationException("Pipeline ${name} has no published release")
        if( pipeline.revoked || release.status == PipelineReleaseStatus.REVOKED )
            throw new AbortOperationException("Pipeline ${name}@${release.version} is revoked")

        workflowFile(dir, name, pipeline.workflowPath ?: DEFAULT_WORKFLOW)

        final existing = Files.exists(dir) ? PipelineInfo.load(dir) : null
        final modified = Files.exists(dir) && (existing == null || existing.isModified(dir))
        if( !modified && existing?.version == release.version )
            return existing
        if( modified && !force )
            throw new AbortOperationException("Pipeline ${name} has local modifications -- use -force to overwrite")

        final info = new PipelineInfo(
            name: name,
            version: release.version,
            workflowPath: pipeline.workflowPath ?: DEFAULT_WORKFLOW,
            source: 'git',
            repositoryUrl: pipeline.repositoryUrl,
            revision: release.gitRevision,
            commitId: release.commitId,
            registryUrl: client.primaryUrl())
        fetchInto(dir, info)
        log.info "Installed pipeline ${name}@${release.version}"
        return info
    }

    List<PipelineInfo> list() {
        final root = baseDir.resolve(RemoteModuleResolver.PIPELINES_DIR)
        final result = new ArrayList<PipelineInfo>()
        if( !Files.isDirectory(root) )
            return result
        for( final scope : children(root) ) {
            for( final dir : children(scope) ) {
                final info = PipelineInfo.load(dir)
                if( info != null )
                    result.add(info)
            }
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
        dir.toFile().deleteDir()
        return true
    }

    private void fetchInto(Path dir, PipelineInfo info) {
        Files.createDirectories(dir.parent)
        final tmp = Files.createTempDirectory(dir.parent, ".tmp-${dir.fileName}-")
        try {
            fetcher.fetch(info.repositoryUrl, info.revision, info.commitId, tmp)
            info.checksum = PipelineInfo.checksum(tmp)
            info.save(tmp)
            if( Files.exists(dir) )
                dir.toFile().deleteDir()
            Files.move(tmp, dir, StandardCopyOption.ATOMIC_MOVE)
        }
        catch( Throwable e ) {
            tmp.toFile().deleteDir()
            throw e
        }
    }

    // skips in-flight '.tmp-' dirs
    private static List<Path> children(Path dir) {
        final result = new ArrayList<Path>()
        Files.newDirectoryStream(dir).withCloseable { stream ->
            for( final p : stream ) {
                if( Files.isDirectory(p) && !p.fileName.toString().startsWith('.') )
                    result.add(p)
            }
        }
        return result
    }
}
