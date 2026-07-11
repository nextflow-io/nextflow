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

package nextflow.install2r

import java.nio.file.FileSystems
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
 * Handle install2.r R library creation and caching
 */
@Slf4j
@CompileStatic
class Install2rCache {

    /**
     * Cache the prefix path for each install2.r R library
     */
    static final private Map<String,DataflowVariable<Path>> install2rPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The install2.r settings defined in the nextflow config file
     */
    private Install2rConfig config

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout

    private String installOptions

    private String repos

    private Path configCacheDir0

    @PackageScope String getRepos() { repos }

    @PackageScope String getInstallOptions() { installOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { SysEnv.get() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @TestOnly
    protected Install2rCache() {}

    /**
     * Create an install2.r library cache object
     *
     * @param config An {@link Install2rConfig} object
     */
    Install2rCache(Install2rConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.installOptions() )
            installOptions = config.installOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()

        repos = config.repos() ?: 'https://cloud.r-project.org'
    }

    /**
     * Retrieve the directory where store the install2.r R library.
     *
     * It tries these settings in the following order:
     * 1) {@code install2r.cacheDir} setting in the nextflow config file;
     * 2) the {@code NXF_INSTALL2R_CACHEDIR} environment variable;
     * 3) the {@code $workDir/install2r} path
     *
     * @return
     *      the {@code Path} where store the install2.r libraries
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_INSTALL2R_CACHEDIR )
            cacheDir = getEnv().NXF_INSTALL2R_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('install2r')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store install2.r libraries to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `install2r.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create install2.r cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    /**
     * Get the path on the file system where store an install2.r R library
     *
     * @param spec The install2.r package specification
     * @return the install2.r unique prefix {@link Path} where the library is created
     */
    @PackageScope
    Path install2rPrefixPath(String spec, String installOptionsOverride = null) {
        assert spec

        // it's interpreted as a user provided library directory
        if( spec.contains('/') && (spec as Path).isDirectory() ) {
            return spec as Path
        }

        String content = spec
        // a per-process install-options override yields a distinct environment
        if( installOptionsOverride ) content += "\nopts:$installOptionsOverride"

        final hash = CacheHelper.hasher(content).hash().toString()
        return getCacheDir().resolve("env-$hash")
    }

    /**
     * Run the install2.r tool to create a library and install packages.
     *
     * @param spec The install2.r package definition
     * @param prefixPath The target path for the library
     * @return the install2.r library prefix {@link Path}
     */
    @PackageScope
    Path createLocalInstall2rEnv(String spec, Path prefixPath, String installOptionsOverride = null) {

        if( prefixPath.isDirectory() ) {
            log.debug "install2.r found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the install2.r library $spec -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalInstall2rEnv0(spec, prefixPath, installOptionsOverride) }
        }
        finally {
            file.delete()
        }

        return prefixPath
    }

    @PackageScope
    Path createLocalInstall2rEnv0(String spec, Path prefixPath, String installOptionsOverride = null) {
        if( prefixPath.isDirectory() ) {
            log.debug "install2.r found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        log.info "Creating env using install2.r: $spec [cache $prefixPath]"

        // per-process `options` override the config-level `install2r.installOptions`
        final effectiveOptions = installOptionsOverride ?: installOptions
        def opts = effectiveOptions ? "${effectiveOptions} " : ''
        def repoOpt = repos ? "-r ${Escape.cli(repos)} " : ''
        // --error makes a failed package install exit non-zero (install2.r warns
        // but exits 0 by default), so a broken env is detected instead of cached
        def cmd = "mkdir -p ${Escape.path(prefixPath)} && install2.r --error ${repoOpt}-l ${Escape.path(prefixPath)} ${opts}${spec}"

        try {
            runCommand( cmd )
            log.debug "'install2.r' create complete env=$spec path=$prefixPath"
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
        log.trace """install2.r create
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
            def msg = "Failed to create install2.r library\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given an install2.r package string returns a {@link DataflowVariable} which holds
     * the local library path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * library creation is actually executed.
     *
     * @param spec
     *      install2.r package string
     * @return
     *      The {@link DataflowVariable} which hold (and create) the local library
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String spec, String installOptionsOverride = null) {
        final prefixPath = install2rPrefixPath(spec, installOptionsOverride)
        final install2rEnvPath = prefixPath.toString()
        if( install2rEnvPath in install2rPrefixPaths ) {
            log.trace "install2.r found local environment `$spec`"
            return install2rPrefixPaths[install2rEnvPath]
        }

        synchronized (install2rPrefixPaths) {
            def result = install2rPrefixPaths[install2rEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalInstall2rEnv(spec, prefixPath, installOptionsOverride) })
                install2rPrefixPaths[install2rEnvPath] = result
            }
            else {
                log.trace "install2.r found local cache for environment `$spec` (2)"
            }
            return result
        }
    }

    /**
     * Create an install2.r library caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * library is actually created.
     *
     * @param spec The install2.r package string
     * @return the local library path prefix {@link Path}
     */
    Path getCachePathFor(String spec, String installOptionsOverride = null) {
        def promise = getLazyImagePath(spec, installOptionsOverride)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create install2.r library `$spec`")
        log.trace "install2.r cache for env `$spec` path=$result"
        return result
    }

}
