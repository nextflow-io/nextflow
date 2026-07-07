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

package nextflow.nix

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
 * Handle Nix environment creation and caching
 */
@Slf4j
@CompileStatic
class NixCache {

    /**
     * Cache the prefix path for each Nix environment
     */
    static final private Map<String,DataflowVariable<Path>> nixPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The Nix settings defined in the nextflow config file
     */
    private NixConfig config

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout

    private String installOptions

    private String flakeRef

    private Path configCacheDir0

    @PackageScope String getInstallOptions() { installOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { SysEnv.get() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @PackageScope String getFlakeRef() { flakeRef }

    @TestOnly
    protected NixCache() {}

    /**
     * Create a Nix env cache object
     *
     * @param config A {@link NixConfig} object
     */
    NixCache(NixConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.installOptions() )
            installOptions = config.installOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()

        flakeRef = config.flakeRef() ?: 'nixpkgs'
    }

    /**
     * Retrieve the directory where store the Nix environment.
     *
     * It tries these settings in the following order:
     * 1) {@code nix.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/nix} path
     *
     * @return
     *      the {@code Path} where store the Nix envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_NIX_CACHEDIR )
            cacheDir = getEnv().NXF_NIX_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('nix')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Nix environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `nix.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create Nix cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    /**
     * Get the path on the file system where store a Nix environment
     *
     * @param spec The Nix environment specification
     * @return the Nix unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path nixPrefixPath(String spec, String installOptionsOverride = null) {
        assert spec

        // it's interpreted as a user provided profile directory
        if( spec.contains('/') && (spec as Path).isDirectory() ) {
            return spec as Path
        }

        String content = spec + "\nflake:" + flakeRef

        // a per-process install-options override yields a distinct environment
        if( installOptionsOverride ) content += "\nopts:$installOptionsOverride"

        final hash = CacheHelper.hasher(content).hash().toString()
        return getCacheDir().resolve("env-$hash")
    }

    /**
     * Run the Nix tool to create an environment and install packages.
     *
     * @param spec The Nix environment definition
     * @param prefixPath The target path for the environment
     * @return the Nix environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalNixEnv(String spec, Path prefixPath, String installOptionsOverride = null) {

        if( prefixPath.isDirectory() ) {
            log.debug "Nix found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the Nix environment $spec -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalNixEnv0(spec, prefixPath, installOptionsOverride) }
        }
        finally {
            file.delete()
        }

        return prefixPath
    }

    @PackageScope
    Path createLocalNixEnv0(String spec, Path prefixPath, String installOptionsOverride = null) {
        if( prefixPath.isDirectory() ) {
            log.debug "Nix found local env for environment=$spec; path=$prefixPath"
            return prefixPath
        }

        log.info "Creating env using Nix: $spec [cache $prefixPath]"

        // map each token to a flake installable, e.g. 'bwa' -> 'nixpkgs#bwa'
        final installables = spec.tokenize().collect { String t ->
            t.contains('#') ? t : "${flakeRef}#${t}".toString()
        }.join(' ')

        // per-process `options` override the config-level `nix.installOptions`
        final effectiveOptions = installOptionsOverride ?: installOptions
        String opts = effectiveOptions ? "$effectiveOptions " : ''
        def cmd = "nix --extra-experimental-features 'nix-command flakes' profile install --profile ${Escape.path(prefixPath)} ${opts}${installables}"

        try {
            runCommand( cmd )
            log.debug "'nix' create complete env=$spec path=$prefixPath"
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
        log.trace """nix create
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
            def msg = "Failed to create Nix environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given a Nix environment string returns a {@link DataflowVariable} which holds
     * the local environment path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment creation is actually executed.
     *
     * @param spec
     *      Nix environment string
     * @return
     *      The {@link DataflowVariable} which hold (and create) the local environment
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String spec, String installOptionsOverride = null) {
        final prefixPath = nixPrefixPath(spec, installOptionsOverride)
        final nixEnvPath = prefixPath.toString()
        if( nixEnvPath in nixPrefixPaths ) {
            log.trace "Nix found local environment `$spec`"
            return nixPrefixPaths[nixEnvPath]
        }

        synchronized (nixPrefixPaths) {
            def result = nixPrefixPaths[nixEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalNixEnv(spec, prefixPath, installOptionsOverride) })
                nixPrefixPaths[nixEnvPath] = result
            }
            else {
                log.trace "Nix found local cache for environment `$spec` (2)"
            }
            return result
        }
    }

    /**
     * Create a Nix environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param spec The Nix environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String spec, String installOptionsOverride = null) {
        def promise = getLazyImagePath(spec, installOptionsOverride)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create Nix environment `$spec`")
        log.trace "Nix cache for env `$spec` path=$result"
        return result
    }

}
