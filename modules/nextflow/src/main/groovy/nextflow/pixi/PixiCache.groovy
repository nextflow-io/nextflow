/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.pixi

import java.nio.file.FileSystems
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.LazyDataflowVariable
import nextflow.Global
import nextflow.file.FileMutex
import nextflow.util.CacheHelper
import nextflow.util.Duration
import nextflow.util.Escape
import nextflow.util.TestOnly

/**
 * Handle Pixi environment creation and caching
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
@Slf4j
@CompileStatic
class PixiCache {

    /**
     * Cache the prefix path for each Pixi environment
     */
    static final private Map<String,DataflowVariable<Path>> pixiPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The Pixi settings defined in the nextflow config file
     */
    private PixiConfig config

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout = Duration.of('20min')

    private String createOptions

    private Path configCacheDir0

    @PackageScope String getCreateOptions() { createOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { System.getenv() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @TestOnly
    protected PixiCache() {}

    /**
     * Create a Pixi env cache object
     *
     * @param config A {@link PixiConfig} object
     */
    PixiCache(PixiConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.createOptions() )
            createOptions = config.createOptions()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir().toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the pixi environment.
     *
     * If tries these setting in the following order:
     * 1) {@code pixi.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/pixi} path
     *
     * @return
     *      the {@code Path} where store the pixi envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_PIXI_CACHEDIR )
            cacheDir = getEnv().NXF_PIXI_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('pixi')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Pixi environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `pixi.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create Pixi cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    @PackageScope
    boolean isTomlFilePath(String str) {
        str.endsWith('.toml') && !str.contains('\n')
    }

    @PackageScope
    boolean isLockFilePath(String str) {
        str.endsWith('.lock') && !str.contains('\n')
    }

    /**
     * Get the path on the file system where store a Pixi environment
     *
     * @param pixiEnv The pixi environment
     * @return the pixi unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path pixiPrefixPath(String pixiEnv) {
        assert pixiEnv

        String content
        String name = 'env'

        // check if it's a TOML file (pixi.toml or pyproject.toml)
        if( isTomlFilePath(pixiEnv) ) {
            try {
                final path = pixiEnv as Path
                content = path.text
                name = path.baseName
            }
            catch( NoSuchFileException e ) {
                throw new IllegalArgumentException("Pixi environment file does not exist: $pixiEnv")
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error parsing Pixi environment TOML file: $pixiEnv -- Check the log file for details", e)
            }
        }
        // check if it's a lock file (pixi.lock)
        else if( isLockFilePath(pixiEnv) ) {
            try {
                final path = pixiEnv as Path
                content = path.text
                name = path.baseName
            }
            catch( NoSuchFileException e ) {
                throw new IllegalArgumentException("Pixi lock file does not exist: $pixiEnv")
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error parsing Pixi lock file: $pixiEnv -- Check the log file for details", e)
            }
        }
        // it's interpreted as user provided prefix directory
        else if( pixiEnv.contains('/') ) {
            final prefix = pixiEnv as Path
            if( !prefix.isDirectory() )
                throw new IllegalArgumentException("Pixi prefix path does not exist or is not a directory: $prefix")
            if( prefix.fileSystem != FileSystems.default )
                throw new IllegalArgumentException("Pixi prefix path must be a POSIX file path: $prefix")

            return prefix
        }
        else if( pixiEnv.contains('\n') ) {
            throw new IllegalArgumentException("Invalid Pixi environment definition: $pixiEnv")
        }
        else {
            // it's interpreted as a package specification
            content = pixiEnv
        }

        final hash = CacheHelper.hasher(content).hash().toString()
        getCacheDir().resolve("$name-$hash")
    }

    /**
     * Run the pixi tool to create an environment in the file system.
     *
     * @param pixiEnv The pixi environment definition
     * @return the pixi environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalPixiEnv(String pixiEnv, Path prefixPath) {

        if( prefixPath.isDirectory() ) {
            log.debug "pixi found local env for environment=$pixiEnv; path=$prefixPath"
            return prefixPath
        }

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the pixi environment $pixiEnv -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalPixiEnv0(pixiEnv, prefixPath) }
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
    Path createLocalPixiEnv0(String pixiEnv, Path prefixPath) {
        log.info "Creating env using pixi: $pixiEnv [cache $prefixPath]"

        String opts = createOptions ? "$createOptions " : ''

        def cmd
        if( isTomlFilePath(pixiEnv) || isLockFilePath(pixiEnv) ) {
            final target = Escape.path(makeAbsolute(pixiEnv))
            final projectDir = makeAbsolute(pixiEnv).parent

            // Create environment from project file
            cmd = "cd ${Escape.path(projectDir)} && pixi install ${opts}"

            // Set up the environment directory
            prefixPath.mkdirs()
            final envLink = prefixPath.resolve('.pixi')
            if( !envLink.exists() ) {
                envLink.toFile().createNewFile()
                envLink.write(projectDir.toString())
            }
        }
        else {
            // Create environment from package specification
            prefixPath.mkdirs()
            final manifestFile = prefixPath.resolve('pixi.toml')

            // Create a simple pixi.toml with the requested packages
            manifestFile.text = """\
[project]
name = "nextflow-env"
version = "0.1.0"
description = "Nextflow generated Pixi environment"
channels = ["conda-forge"]

[dependencies]
${pixiEnv}
""".stripIndent()

            cmd = "cd ${Escape.path(prefixPath)} && pixi install ${opts}"
        }

        try {
            runCommand( cmd )
            log.debug "'pixi' create complete env=$pixiEnv path=$prefixPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted image file
            prefixPath.delete()
            throw e
        }
        return prefixPath
    }

    @PackageScope
    int runCommand( String cmd ) {
        log.trace """pixi create
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
            def msg = "Failed to create Pixi environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    /**
     * Given a remote image URL returns a {@link DataflowVariable} which holds
     * the local image path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * image download is actually executed.
     *
     * @param pixiEnv
     *      Pixi environment string
     * @return
     *      The {@link DataflowVariable} which hold (and pull) the local image file
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String pixiEnv) {
        final prefixPath = pixiPrefixPath(pixiEnv)
        final pixiEnvPath = prefixPath.toString()
        if( pixiEnvPath in pixiPrefixPaths ) {
            log.trace "pixi found local environment `$pixiEnv`"
            return pixiPrefixPaths[pixiEnvPath]
        }

        synchronized (pixiPrefixPaths) {
            def result = pixiPrefixPaths[pixiEnvPath]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalPixiEnv(pixiEnv, prefixPath) })
                pixiPrefixPaths[pixiEnvPath] = result
            }
            else {
                log.trace "pixi found local cache for environment `$pixiEnv` (2)"
            }
            return result
        }
    }

    /**
     * Create a pixi environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param pixiEnv The pixi environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String pixiEnv) {
        def promise = getLazyImagePath(pixiEnv)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create Pixi environment `$pixiEnv`")
        log.trace "Pixi cache for env `$pixiEnv` path=$result"
        return result
    }

}
