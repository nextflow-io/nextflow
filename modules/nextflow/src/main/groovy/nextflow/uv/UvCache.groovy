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

package nextflow.uv

import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.Paths
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
 * Handle uv virtual environment creation and caching
 *
 * @author Evan Floden
 */
@Slf4j
@CompileStatic
class UvCache {

    /**
     * Cache the prefix path for each uv environment
     */
    static final private Map<String,DataflowVariable<Path>> uvPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The uv settings defined in the nextflow config file
     */
    private UvConfig config

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout

    private String installOptions

    private String pythonVersion

    private Path configCacheDir0

    @PackageScope String getInstallOptions() { installOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { SysEnv.get() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @PackageScope String getPythonVersion() { pythonVersion }

    @TestOnly
    protected UvCache() {}

    /**
     * Create a uv env cache object
     *
     * @param config A {@link UvConfig} object
     */
    UvCache(UvConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.installOptions() )
            installOptions = config.installOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()

        if( config.pythonVersion() )
            pythonVersion = config.pythonVersion()
    }

    /**
     * Retrieve the directory where store the uv environment.
     *
     * It tries these settings in the following order:
     * 1) {@code uv.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/uv} path
     *
     * @return
     *      the {@code Path} where store the uv envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_UV_CACHEDIR )
            cacheDir = getEnv().NXF_UV_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('uv')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store uv environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `uv.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create uv cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    @PackageScope
    boolean isRequirementsFile(String str) {
        (str.endsWith('.txt') || str.endsWith('.in')) && !str.contains('\n')
    }

    @PackageScope
    boolean isPyProjectFile(String str) {
        str.endsWith('pyproject.toml') && !str.contains('\n')
    }

    /**
     * Get the path on the file system where store a uv environment
     *
     * @param uvEnv The uv environment specification
     * @return the uv unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path uvPrefixPath(String uvEnv) {
        assert uvEnv

        String content
        String name = 'env'
        // check if it's a requirements file
        if( isRequirementsFile(uvEnv) || isPyProjectFile(uvEnv) ) {
            try {
                final path = uvEnv as Path
                content = path.text
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error reading uv environment file: $uvEnv -- Check the log file for details", e)
            }
        }
        // it's interpreted as user provided prefix directory
        else if( uvEnv.contains('/') ) {
            final prefix = uvEnv as Path
            if( !prefix.isDirectory() )
                throw new IllegalArgumentException("uv environment path does not exist or is not a directory: $prefix")
            if( prefix.fileSystem != FileSystems.default )
                throw new IllegalArgumentException("uv environment path must be a POSIX file path: $prefix")

            return prefix
        }
        else if( uvEnv.contains('\n') ) {
            throw new IllegalArgumentException("Invalid uv environment definition: $uvEnv")
        }
        else {
            content = uvEnv
        }

        // include python version in hash if specified
        if( pythonVersion ) content += "\npython:$pythonVersion"

        final hash = CacheHelper.hasher(content).hash().toString()
        getCacheDir().resolve("$name-$hash")
    }

    /**
     * Run the uv tool to create a virtual environment and install packages.
     *
     * @param uvEnv The uv environment definition
     * @param prefixPath The target path for the virtual environment
     * @return the uv environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalUvEnv(String uvEnv, Path prefixPath) {

        if( prefixPath.isDirectory() ) {
            log.debug "uv found local env for environment=$uvEnv; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the uv environment $uvEnv -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalUvEnv0(uvEnv, prefixPath) }
        }
        finally {
            file.delete()
        }

        return prefixPath
    }

    @PackageScope
    Path makeAbsolute( String envFile ) {
        Paths.get(envFile).toAbsolutePath()
    }

    @PackageScope
    Path createLocalUvEnv0(String uvEnv, Path prefixPath) {
        if( prefixPath.isDirectory() ) {
            log.debug "uv found local env for environment=$uvEnv; path=$prefixPath"
            return prefixPath
        }

        log.info "Creating env using uv: $uvEnv [cache $prefixPath]"

        String opts = installOptions ? "$installOptions " : ''
        String pythonOpt = pythonVersion ? "--python $pythonVersion " : ''

        def cmd
        // First create the virtual environment
        cmd = "uv venv ${pythonOpt}${Escape.path(prefixPath)} && "

        if( isRequirementsFile(uvEnv) ) {
            cmd += "uv pip install ${opts}--python ${Escape.path(prefixPath.resolve('bin/python'))} -r ${Escape.path(makeAbsolute(uvEnv))}"
        }
        else if( isPyProjectFile(uvEnv) ) {
            cmd += "uv pip install ${opts}--python ${Escape.path(prefixPath.resolve('bin/python'))} -r ${Escape.path(makeAbsolute(uvEnv))}"
        }
        else {
            // space-separated package list e.g. 'numpy pandas matplotlib'
            cmd += "uv pip install ${opts}--python ${Escape.path(prefixPath.resolve('bin/python'))} $uvEnv"
        }

        try {
            runCommand( cmd )
            log.debug "'uv' create complete env=$uvEnv path=$prefixPath"
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
        log.trace """uv create
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
            def msg = "Failed to create uv environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given a uv environment string returns a {@link DataflowVariable} which holds
     * the local environment path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment creation is actually executed.
     *
     * @param uvEnv
     *      uv environment string
     * @return
     *      The {@link DataflowVariable} which hold (and create) the local environment
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String uvEnv) {
        final prefixPath = uvPrefixPath(uvEnv)
        final uvEnvPath = prefixPath.toString()
        if( uvEnvPath in uvPrefixPaths ) {
            log.trace "uv found local environment `$uvEnv`"
            return uvPrefixPaths[uvEnvPath]
        }

        synchronized (uvPrefixPaths) {
            def result = uvPrefixPaths[uvEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalUvEnv(uvEnv, prefixPath) })
                uvPrefixPaths[uvEnvPath] = result
            }
            else {
                log.trace "uv found local cache for environment `$uvEnv` (2)"
            }
            return result
        }
    }

    /**
     * Create a uv environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param uvEnv The uv environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String uvEnv) {
        def promise = getLazyImagePath(uvEnv)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create uv environment `$uvEnv`")
        log.trace "uv cache for env `$uvEnv` path=$result"
        return result
    }

}
