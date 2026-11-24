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
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.concurrent.ConcurrentHashMap
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.Const
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.module.spi.RemoteModuleResolver

/**
 * Resolves, installs, lists and removes registry pipelines stored under {@code <baseDir>/pipelines}
 */
@Slf4j
@CompileStatic
class PipelineManager {

    private static final Pattern NAME_PATTERN = ~/^[a-z0-9][a-z0-9._-]*\/[a-z][a-z0-9._-]*$/

    /** An installed pipeline directory with its {@code .pipeline-info} */
    static record InstalledPipeline(Path dir, PipelineInfo info) {}

    /** Outcome of {@link #install}: {@code alreadyInstalled} means nothing was fetched */
    static record InstallResult(PipelineInfo info, boolean alreadyInstalled) {}

    // never evicted: fine for a CLI or run lifetime, a long-lived process would need a TTL
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
        if( info != null ) {
            if( isModified(info, dir) )
                log.warn1 "Pipeline ${name} has local modifications (checksum mismatch)"
        }
        else if( Files.exists(dir) )
            log.warn1 "Pipeline ${name} has no registry origin (.pipeline-info missing)"
        else
            install(name, null, false)
        return dir.resolve(Const.DEFAULT_MAIN_FILE_NAME)
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

    InstallResult install(String name, String version, boolean force) {
        final dir = pipelineDir(name)
        final info = resolveRelease(client, name, version)

        final exists = Files.exists(dir)
        final existing = exists ? PipelineInfo.load(dir) : null
        final modified = exists && (existing == null || existing.isModified(dir))
        if( !modified && existing?.version == info.version )
            return new InstallResult(existing, true)
        if( modified && !force )
            throw new AbortOperationException("Pipeline ${name} has local modifications -- use -force to overwrite")

        log.info "Installing pipeline ${name}@${info.version}..."
        fetchInto(dir, info)
        log.info "Pipeline ${name}@${info.version} installed in ${dir.toAbsolutePath()}"
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
            final commit = fetcher.fetch(info.repositoryUrl, info.revision, info.commitId, tmp)
            if( !Files.isRegularFile(tmp.resolve(Const.DEFAULT_MAIN_FILE_NAME)) )
                throw new AbortOperationException("Pipeline ${info.name}@${info.version} has no main.nf -- a registry pipeline must have a top-level main.nf")
            info.commitId = info.commitId ?: commit
            info.checksum = PipelineInfo.checksum(tmp)
            info.save(tmp)
            MODIFIED.remove(realPath(dir))
            replace(dir, tmp)
        }
        catch( Throwable e ) {
            try {
                FileHelper.deletePath(tmp)
            }
            catch( IOException ex ) {
                log.warn "Cannot delete temporary directory ${tmp}: ${ex.message}"
            }
            throw e
        }
    }

    // the existing copy is moved aside, and moved back when the new one can't be moved in
    private void replace(Path dir, Path tmp) {
        final old = Files.exists(dir, LinkOption.NOFOLLOW_LINKS) ? dir.resolveSibling("${tmp.fileName}-old") : null
        if( old )
            move(dir, old)
        try {
            move(tmp, dir)
        }
        catch( Throwable e ) {
            if( old )
                move(old, dir)
            throw e
        }
        if( old )
            FileHelper.deletePath(old)
    }

    protected void move(Path source, Path target) {
        Files.move(source, target, StandardCopyOption.ATOMIC_MOVE)
    }

    /** Thrown when the registry has no pipeline of the given name */
    @InheritConstructors
    static class PipelineNotFoundException extends AbortOperationException {}

    /** A registry pipeline and one of its releases */
    static record RegistryRelease(Pipeline pipeline, PipelineRelease release) {}

    /**
     * Pipeline and release from the registry, the latest release when no version is given.
     * Fails with {@link PipelineNotFoundException} when the pipeline is missing, and with {@link AbortOperationException} on any other error.
     */
    static RegistryRelease fetchRelease(RegistryClient client, String name, String version) {
        if( !isValidName(name) )
            throw new AbortOperationException("Invalid pipeline name '$name' -- expected <scope>/<name>, e.g. nf-core/rnaseq")
        Pipeline pipeline = null
        try {
            pipeline = client.getPipeline(name)
        }
        catch( RegistryException e ) {
            if( e.statusCode != 404 )
                throw registryError(client, e)
        }
        if( pipeline == null )
            throw new PipelineNotFoundException("Pipeline ${name} not found in the registry")
        PipelineRelease release = pipeline.latest
        if( version ) {
            try {
                release = client.getPipelineRelease(name, version)
            }
            catch( RegistryException e ) {
                if( e.statusCode != 404 )
                    throw registryError(client, e)
                release = null
            }
            if( release == null )
                throw new AbortOperationException("Version ${version} of pipeline ${name} not found in the registry")
        }
        if( release == null )
            throw new AbortOperationException("Pipeline ${name} has no published release")
        return new RegistryRelease(pipeline, release)
    }

    /** A registry error other than not found */
    static AbortOperationException registryError(RegistryClient client, RegistryException e) {
        return new AbortOperationException("Cannot query registry ${client.primaryUrl()}: ${e.message}", e)
    }

    /** Release of a registry pipeline, the latest when no version is given; fails when it is missing or revoked */
    static PipelineInfo resolveRelease(RegistryClient client, String name, String version) {
        final result = fetchRelease(client, name, version)
        final pipeline = result.pipeline()
        final release = result.release()
        if( pipeline.revoked || release.status == PipelineReleaseStatus.REVOKED )
            throw new AbortOperationException("Pipeline ${name}@${release.version} is revoked")
        return new PipelineInfo(
            name: name,
            version: release.version,
            repositoryUrl: pipeline.repositoryUrl,
            revision: release.gitRevision,
            commitId: release.commitId,
            registryUrl: client.primaryUrl())
    }
}
