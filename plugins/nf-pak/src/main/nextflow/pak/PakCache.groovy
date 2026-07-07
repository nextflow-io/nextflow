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

package nextflow.pak

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
 * Handle R pak environment creation and caching
 */
@Slf4j
@CompileStatic
class PakCache {

    /**
     * Cache the prefix path for each R pak environment
     */
    static final private Map<String,DataflowVariable<Path>> pakPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The pak settings defined in the nextflow config file
     */
    private PakConfig config

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
    protected PakCache() {}

    /**
     * Create a pak env cache object
     *
     * @param config A {@link PakConfig} object
     */
    PakCache(PakConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.installOptions() )
            installOptions = config.installOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the R pak environment.
     *
     * It tries these settings in the following order:
     * 1) {@code pak.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/pak} path
     *
     * @return
     *      the {@code Path} where store the pak envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_PAK_CACHEDIR )
            cacheDir = getEnv().NXF_PAK_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('pak')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store R pak environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `pak.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create R pak cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    /**
     * Get the path on the file system where store an R pak environment
     *
     * @param spec The pak environment specification
     * @return the pak unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path pakPrefixPath(String spec, String installOptionsOverride = null) {
        assert spec

        // it's interpreted as a user provided library directory
        if( spec.contains('/') && (spec as Path).isDirectory() ) {
            return spec as Path
        }

        // for a manifest file (renv.lock / DESCRIPTION) hash its content so the
        // cache is invalidated when the file changes
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
     * Run the pak tool to create an environment and install packages.
     *
     * @param spec The pak environment definition
     * @param prefixPath The target path for the environment
     * @return the pak environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalPakEnv(String spec, Path prefixPath, String installOptionsOverride = null) {

        if( prefixPath.isDirectory() ) {
            log.debug "pak found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the R pak environment $spec -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalPakEnv0(spec, prefixPath, installOptionsOverride) }
        }
        finally {
            file.delete()
        }

        return prefixPath
    }

    @PackageScope
    Path createLocalPakEnv0(String spec, Path prefixPath, String installOptionsOverride = null) {
        if( prefixPath.isDirectory() ) {
            log.debug "pak found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        log.info "Creating env using R pak: $spec [cache $prefixPath]"

        // per-process `options` override the config-level `pak.installOptions`
        final effectiveOptions = installOptionsOverride ?: installOptions
        def opts = effectiveOptions ? ", ${effectiveOptions}" : ''
        def cmd
        if( isManifestFile(spec) ) {
            // a DESCRIPTION file: install the declared dependencies of the package
            // in its directory (ask=FALSE for non-interactive safety)
            final path = spec as Path
            cmd = "mkdir -p ${Escape.path(prefixPath)} && Rscript -e 'pak::local_install_deps(\"${path.parent}\", lib=\"${prefixPath}\", ask = FALSE${opts})'"
        }
        else {
            // build an R character vector of package names, e.g. "dplyr", "ggplot2"
            def pkgs = spec.tokenize().collect { "\"" + it + "\"" }.join(', ')
            cmd = "mkdir -p ${Escape.path(prefixPath)} && Rscript -e 'pak::pkg_install(c(${pkgs}), lib=\"${prefixPath}\"${opts})'"
        }

        try {
            runCommand( cmd )
            log.debug "'pak' create complete env=$spec path=$prefixPath"
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
        log.trace """pak create
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
            def msg = "Failed to create R pak environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given a pak environment string returns a {@link DataflowVariable} which holds
     * the local environment path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment creation is actually executed.
     *
     * @param spec
     *      pak environment string
     * @return
     *      The {@link DataflowVariable} which hold (and create) the local environment
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String spec, String installOptionsOverride = null) {
        final prefixPath = pakPrefixPath(spec, installOptionsOverride)
        final pakEnvPath = prefixPath.toString()
        if( pakEnvPath in pakPrefixPaths ) {
            log.trace "pak found local environment `$spec`"
            return pakPrefixPaths[pakEnvPath]
        }

        synchronized (pakPrefixPaths) {
            def result = pakPrefixPaths[pakEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalPakEnv(spec, prefixPath, installOptionsOverride) })
                pakPrefixPaths[pakEnvPath] = result
            }
            else {
                log.trace "pak found local cache for environment `$spec` (2)"
            }
            return result
        }
    }

    /**
     * Create a pak environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param spec The pak environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String spec, String installOptionsOverride = null) {
        def promise = getLazyImagePath(spec, installOptionsOverride)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create R pak environment `$spec`")
        log.trace "pak cache for env `$spec` path=$result"
        return result
    }

}
