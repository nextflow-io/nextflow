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
import org.yaml.snakeyaml.Yaml

/**
 * Handle Pixi environment creation and caching
 *
 * @author GitHub Copilot
 */
@Slf4j
@CompileStatic
class PixiCache {
    static final private Object pixiLock = new Object()

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

    private List<String> channels = Collections.emptyList()

    @PackageScope String getCreateOptions() { createOptions }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { System.getenv() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    @PackageScope List<String> getChannels() { channels }

    @TestOnly
    protected PixiCache() {}

    /**
     * Create a Pixi env cache object
     *
     * @param config A {@link Map} object
     */
    PixiCache(PixiConfig config) {
        this.config = config

        if( config.createTimeout() )
            createTimeout = config.createTimeout()

        if( config.createOptions() )
            createOptions = config.createOptions()

        if( config.getChannels() )
            channels = config.getChannels()

        if( config.cacheDir() )
            configCacheDir0 = config.cacheDir()
    }

    /**
     * Retrieve the directory where to store the pixi environment.
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
        // check if it's a TOML file
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
    Path createLocalPixiEnv(String pixiEnv) {
        final prefixPath = pixiPrefixPath(pixiEnv)

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the pixi environment $pixiEnv -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)

        if( prefixPath.isDirectory() ) {
            log.debug "pixi found local env for environment=$pixiEnv; path=$prefixPath"
            try {
                // have to check environment packages, because they are stored externally in the host
                // and might get modified there
                mutex.lock { checkLocalPixiEnv0(pixiEnv, prefixPath) }
            }
            finally {
                file.delete()
            }
            return prefixPath
        }

        try {
            mutex.lock { createLocalPixiEnv0(pixiEnv, prefixPath) }
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

    Path checkLocalPixiEnv0(String pixiEnv, Path prefixPath) {
        log.debug "Checking env using pixi: $pixiEnv [cache $prefixPath]"

        def cmd = "cd ${Escape.path(prefixPath)} && pixi install"

        try {
            runCommand( cmd )
            log.debug "'pixi' check complete env=$pixiEnv path=$prefixPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted image file
            prefixPath.delete()
            throw e
        }
        return prefixPath
    }

    @PackageScope
    Path createLocalPixiEnv0(String pixiEnv, Path prefixPath) {

        log.info "Creating env using pixi: $pixiEnv [cache $prefixPath]"

        def cmd
        if( isTomlFilePath(pixiEnv) ) {
            cmd = "mkdir -p ${Escape.path(prefixPath)} && cp ${Escape.path(makeAbsolute(pixiEnv))} ${Escape.path(prefixPath)}/pixi.toml && cd ${Escape.path(prefixPath)} && pixi install"
        }
        else {
            // if not a toml file, assume it's a list of packages
            cmd = "mkdir -p ${Escape.path(prefixPath)} && cd ${Escape.path(prefixPath)} && echo '[project]' > pixi.toml && echo 'name = \"nf-pixi-env\"' >> pixi.toml && echo '' >> pixi.toml && echo 'channels = [\"conda-forge\"]' >> pixi.toml && echo 'platforms = [\"linux-64\", \"osx-arm64\", \"win-64\"]' >> pixi.toml && echo '' >> pixi.toml && echo '[dependencies]' >> pixi.toml"
            
            // add each package to the toml file
            pixiEnv.split(' ').each { pkg ->
                cmd += " && echo '$pkg = \"*\"' >> pixi.toml"
            }
            
            // finally run pixi install
            cmd += " && pixi install"
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
    void runCommand( String cmd ) {
        log.trace """pixi env create
                     command: $cmd
                     timeout: $createTimeout""".stripIndent(true)

        final max = createTimeout.toMillis()
        final builder = new ProcessBuilder(['bash','-c',cmd])
        final proc = builder.start()
        final err = new StringBuilder()
        final consumer = proc.consumeProcessErrorStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            def msg = "Failed to create pixi environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
    }

    /**
     * Given a remote image URL returns a {@link DataflowVariable} which holds
     * the local image path.
     *
     * This method synchronize multiple concurrent requests so that only one
     * image download is actually executed.
     *
     * @param pixiEnv
     *      Pixi environment string
     * @return
     *      The {@link DataflowVariable} which hold (and pull) the local image file
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String pixiEnv) {
        if( pixiEnv in pixiPrefixPaths ) {
            log.trace "pixi found local environment `$pixiEnv`"
            return pixiPrefixPaths[pixiEnv]
        }

        synchronized (pixiPrefixPaths) {
            def result = pixiPrefixPaths[pixiEnv]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalPixiEnv(pixiEnv) })
                pixiPrefixPaths[pixiEnv] = result
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
     * This method synchronize multiple concurrent requests so that only one
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
