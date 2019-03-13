/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.conda

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
import org.yaml.snakeyaml.Yaml
/**
 * Handle Conda environment creation and caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CondaCache {

    /**
     * Cache the prefix path for each Conda environment
     */
    static final private Map<String,DataflowVariable<Path>> condaPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The Conda settings defined in the nextflow config file
     */
    private CondaConfig config

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

    /** Only for debugging purpose - do not use */
    @PackageScope
    CondaCache() {}

    /**
     * Create a Conda env cache object
     *
     * @param config A {@link Map} object
     */
    CondaCache(CondaConfig config) {
        this.config = config

        if( config.createTimeout )
            createTimeout = config.createTimeout as Duration

        if( config.createOptions )
            createOptions = config.createOptions

        if( config.cacheDir )
            configCacheDir0 = (config.cacheDir as Path).toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the conda environment.
     *
     * If tries these setting in the following order:
     * 1) {@code conda.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/conda} path
     *
     * @return
     *      the {@code Path} where store the conda envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_CONDA_CACHEDIR )
            cacheDir = getEnv().NXF_CONDA_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('conda')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Conda environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `conda.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create Conda cache directory: $cacheDir -- Make sure a file with the same does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    @PackageScope
    boolean isYamlFilePath(String str) {
        (str.endsWith('.yml') || str.endsWith('.yaml')) && !str.contains('\n')
    }

    boolean isTextFilePath(String str) {
        str.endsWith('.txt') && !str.contains('\n')
    }


    /**
     * Get the path on the file system where store a Conda environment
     *
     * @param condaEnv The conda environment
     * @return the conda unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path condaPrefixPath(String condaEnv) {
        assert condaEnv

        String content
        String name = 'env'
        // check if it's a YAML file
        if( isYamlFilePath(condaEnv) ) {
            try {
                final path = condaEnv as Path
                content = path.text
                final yaml = (Map)new Yaml().load(content)
                if( yaml.name )
                    name = yaml.name
                else
                    name = path.baseName
            }
            catch( NoSuchFileException e ) {
                throw new IllegalArgumentException("Conda environment file does not exist: $condaEnv")
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error parsing Conda environment YAML file: $condaEnv -- Chech the log file for details", e)
            }
        }
        else if( isTextFilePath(condaEnv) )  {
            try {
                final path = condaEnv as Path
                content = path.text
                name = path.baseName
            }
            catch( NoSuchFileException e ) {
                throw new IllegalArgumentException("Conda environment file does not exist: $condaEnv")
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error parsing Conda environment text file: $condaEnv -- Chech the log file for details", e)
            }
        }
        // it's interpreted as user provided prefix directory
        else if( condaEnv.contains('/') ) {
            final prefix = condaEnv as Path
            if( !prefix.isDirectory() )
                throw new IllegalArgumentException("Conda prefix path does not exist or it's not a directory: $prefix")
            if( prefix.fileSystem != FileSystems.default )
                throw new IllegalArgumentException("Conda prefix path must be a POSIX file path: $prefix")

            return prefix
        }
        else if( condaEnv.contains('\n') ) {
            throw new IllegalArgumentException("Invalid Conda environment definition: $condaEnv")
        }
        else {
            content = condaEnv
        }

        final hash = CacheHelper.hasher(content).hash().toString()
        getCacheDir().resolve("$name-$hash")
    }

    /**
     * Run the conda tool to create an environment in the file system.
     *
     * @param condaEnv The conda environment definition
     * @return the conda environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalCondaEnv(String condaEnv) {
        final prefixPath = condaPrefixPath(condaEnv)
        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creatign the Conda environment $condaEnv -- please wait it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { createLocalCondaEnv0(condaEnv, prefixPath) }
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
    Path createLocalCondaEnv0(String condaEnv, Path prefixPath) {

        if( prefixPath.isDirectory() ) {
            log.debug "Conda found local env for environment=$condaEnv; path=$prefixPath"
            return prefixPath
        }
        log.info "Creating Conda env: $condaEnv [cache $prefixPath]"

        final opts = createOptions ? "$createOptions " : ''
        def cmd
        if( isYamlFilePath(condaEnv) ) {
            cmd = "conda env create --prefix ${Escape.path(prefixPath)} --file ${Escape.path(makeAbsolute(condaEnv))}"
        }
        else if( isTextFilePath(condaEnv) ) {
            cmd = "conda create $opts--mkdir --yes --quiet --prefix ${Escape.path(prefixPath)} --file ${Escape.path(makeAbsolute(condaEnv))}"
        }
        else {
            cmd = "conda create $opts--mkdir --yes --quiet --prefix ${Escape.path(prefixPath)} $condaEnv"
        }

        try {
            runCommand( cmd )
            log.debug "Conda create complete env=$condaEnv path=$prefixPath"
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
        log.trace """Conda create
                     command: $cmd
                     timeout: $createTimeout""".stripIndent()

        final max = createTimeout.toMillis()
        final builder = new ProcessBuilder(['bash','-c',cmd])
        final proc = builder.start()
        final err = new StringBuilder()
        final consumer = proc.consumeProcessErrorStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            def msg = "Failed to create Conda environment\n  command: $cmd\n  status : $status\n  message:\n"
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
     * @param condaEnv
     *      Conda environment string
     * @return
     *      The {@link DataflowVariable} which hold (and pull) the local image file
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String condaEnv) {
        if( condaEnv in condaPrefixPaths ) {
            log.trace "Conda found local environment `$condaEnv`"
            return condaPrefixPaths[condaEnv]
        }

        synchronized (condaPrefixPaths) {
            def result = condaPrefixPaths[condaEnv]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalCondaEnv(condaEnv) })
                condaPrefixPaths[condaEnv] = result
            }
            else {
                log.trace "Conda found local cache for environment `$condaEnv` (2)"
            }
            return result
        }
    }

    /**
     * Create a conda environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param condaEnv The conda environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String condaEnv) {
        def promise = getLazyImagePath(condaEnv)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create Conda environment `$condaEnv`")
        log.trace "Conda cache for env `$condaEnv` path=$result"
        return result
    }

}
