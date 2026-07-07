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

package nextflow.guix

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.LazyDataflowVariable
import nextflow.Global
import nextflow.SysEnv
import nextflow.file.FileMutex
import nextflow.util.CacheHelper
import nextflow.util.Duration
import nextflow.util.Escape
import nextflow.util.TestOnly

/**
 * Handle GNU Guix environment creation and caching
 */
@Slf4j
@CompileStatic
class GuixCache {

    /**
     * Cache the prefix path for each GNU Guix environment
     */
    static final private Map<String,DataflowVariable<Path>> guixPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The GNU Guix settings defined in the nextflow config file
     */
    private GuixConfig config

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout

    private String installOptions

    private Path configCacheDir0

    @PackageScope String getInstallOptions() { installOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { SysEnv.get() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @TestOnly
    protected GuixCache() {}

    /**
     * Create a GNU Guix env cache object
     *
     * @param config A {@link GuixConfig} object
     */
    GuixCache(GuixConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.installOptions() )
            installOptions = config.installOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the GNU Guix environment.
     *
     * It tries these settings in the following order:
     * 1) {@code guix.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/guix} path
     *
     * @return
     *      the {@code Path} where store the GNU Guix envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_GUIX_CACHEDIR )
            cacheDir = getEnv().NXF_GUIX_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('guix')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store GNU Guix environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `guix.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create GNU Guix cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    /**
     * Get the path on the file system where store a GNU Guix environment
     *
     * @param spec The GNU Guix environment specification
     * @return the GNU Guix unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path guixPrefixPath(String spec, String installOptionsOverride = null) {
        assert spec

        // it's interpreted as user provided prefix directory
        if( spec.contains('/') && (spec as Path).isDirectory() ) {
            return spec as Path
        }

        // for a manifest file (manifest.scm) hash its content so the cache is
        // invalidated when the file changes
        String content = isManifestFile(spec) ? (spec as Path).text : spec

        // a per-process install-options override yields a distinct environment
        if( installOptionsOverride ) content += "\nopts:$installOptionsOverride"

        final hash = CacheHelper.hasher(content).hash().toString()
        return getCacheDir().resolve("env-$hash")
    }

    @PackageScope
    boolean isManifestFile(String spec) {
        spec.contains('/') && Files.isRegularFile(spec as Path)
    }

    /**
     * Run the Guix tool to create an environment and install packages.
     *
     * @param spec The GNU Guix environment definition
     * @param prefixPath The target path for the environment
     * @return the GNU Guix environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalGuixEnv(String spec, Path prefixPath, String installOptionsOverride = null) {

        if( prefixPath.isDirectory() ) {
            log.debug "GNU Guix found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the GNU Guix environment $spec -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalGuixEnv0(spec, prefixPath, installOptionsOverride) }
        }
        finally {
            file.delete()
        }

        return prefixPath
    }

    @PackageScope
    Path createLocalGuixEnv0(String spec, Path prefixPath, String installOptionsOverride = null) {
        if( prefixPath.isDirectory() ) {
            log.debug "GNU Guix found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        log.info "Creating env using GNU Guix: $spec [cache $prefixPath]"

        // per-process `options` override the config-level `guix.installOptions`
        final effectiveOptions = installOptionsOverride ?: installOptions
        String opts = effectiveOptions ? "$effectiveOptions " : ''
        def cmd = isManifestFile(spec)
            ? "guix package --profile=${Escape.path(prefixPath)} ${opts}--manifest=${Escape.path(spec as Path)}"
            : "guix package --profile=${Escape.path(prefixPath)} ${opts}--install $spec"

        try {
            runCommand( cmd )
            log.debug "'guix' create complete env=$spec path=$prefixPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted environment
            prefixPath.deleteDir()
            throw e
        }
        return prefixPath
    }

    @PackageScope
    int runCommand( String cmd ) {
        log.trace """guix create
                     command: $cmd
                     timeout: $createTimeout""".stripIndent(true)

        final max = createTimeout.toMillis()
        final builder = new ProcessBuilder(['bash','-c',cmd])
        final proc = builder.redirectErrorStream(true).start()
        final err = new StringBuilder()
        final consumer = proc.consumeProcessOutputStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            def msg = "Failed to create GNU Guix environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given a GNU Guix environment string returns a {@link DataflowVariable} which holds
     * the local environment path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment creation is actually executed.
     *
     * @param spec
     *      GNU Guix environment string
     * @return
     *      The {@link DataflowVariable} which hold (and create) the local environment
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String spec, String installOptionsOverride = null) {
        final prefixPath = guixPrefixPath(spec, installOptionsOverride)
        final guixEnvPath = prefixPath.toString()
        if( guixEnvPath in guixPrefixPaths ) {
            log.trace "GNU Guix found local environment `$spec`"
            return guixPrefixPaths[guixEnvPath]
        }

        synchronized (guixPrefixPaths) {
            def result = guixPrefixPaths[guixEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalGuixEnv(spec, prefixPath, installOptionsOverride) })
                guixPrefixPaths[guixEnvPath] = result
            }
            else {
                log.trace "GNU Guix found local cache for environment `$spec` (2)"
            }
            return result
        }
    }

    /**
     * Create a GNU Guix environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param spec The GNU Guix environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String spec, String installOptionsOverride = null) {
        def promise = getLazyImagePath(spec, installOptionsOverride)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create GNU Guix environment `$spec`")
        log.trace "GNU Guix cache for env `$spec` path=$result"
        return result
    }

}
