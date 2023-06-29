/*
 * Copyright 2022-2023, Pawsey Supercomputing Research Centre
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

package nextflow.spack

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
 * Handle Spack environment creation and caching
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@Slf4j
@CompileStatic
class SpackCache {

    /**
     * Cache the prefix path for each Spack environment
     */
    static final private Map<String,DataflowVariable<Path>> spackPrefixPaths = new ConcurrentHashMap<>()

    /**
     * The Spack settings defined in the nextflow config file
     */
    private SpackConfig config

    private boolean checksum = true

    private Integer parallelBuilds

    /**
     * Timeout after which the environment creation is aborted
     */
    private Duration createTimeout = Duration.of('60min')

    private Path configCacheDir0

    @PackageScope Integer getParallelBuilds() { parallelBuilds }

    @PackageScope Duration getCreateTimeout() { createTimeout }

    @PackageScope Map<String,String> getEnv() { System.getenv() }

    @PackageScope Path getConfigCacheDir0() { configCacheDir0 }

    /** Only for testing purpose - do not use */
    @PackageScope
    SpackCache() {}

    /**
     * Create a Spack env cache object
     *
     * @param config A {@link Map} object
     */
    SpackCache(SpackConfig config) {
        this.config = config

        if( config.checksum )
            checksum = config.checksum as boolean

        if( config.parallelBuilds )
            parallelBuilds = config.parallelBuilds as Integer

        if( config.createTimeout )
            createTimeout = config.createTimeout as Duration

        if( config.cacheDir )
            configCacheDir0 = (config.cacheDir as Path).toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the spack environment.
     *
     * If tries these setting in the following order:
     * 1) {@code spack.cacheDir} setting in the nextflow config file;
     * 2) the {@code $workDir/spack} path
     *
     * @return
     *      the {@code Path} where store the spack envs
     */
    @PackageScope
    Path getCacheDir() {

        def cacheDir = configCacheDir0

        if( !cacheDir && getEnv().NXF_SPACK_CACHEDIR )
            cacheDir = getEnv().NXF_SPACK_CACHEDIR as Path

        if( !cacheDir )
            cacheDir = getSessionWorkDir().resolve('spack')

        if( cacheDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Spack environments to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `spack.cacheDir` config setting")
        }

        if( !cacheDir.exists() && !cacheDir.mkdirs() ) {
            throw new IOException("Failed to create Spack cache directory: $cacheDir -- Make sure a file with the same name does not exist and you have write permission")
        }

        return cacheDir
    }

    @PackageScope Path getSessionWorkDir() {
        Global.session.workDir
    }

    @PackageScope
    boolean isYamlFilePath(String str) {
        str.endsWith('.yaml') && !str.contains('\n')
    }


    /**
     * Get the path on the file system where store a Spack environment
     *
     * @param spackEnv The spack environment
     * @return the spack unique prefix {@link Path} where the env is created
     */
    @PackageScope
    Path spackPrefixPath(String spackEnv, String arch) {
        assert spackEnv

        String content
        String name = 'env'
        // check if it's a YAML file
        if( isYamlFilePath(spackEnv) ) {
            try {
                final path = spackEnv as Path
                content = path.text
                final yaml = (Map)new Yaml().load(content)
                if( yaml.name )
                    name = yaml.name
                else
                    name = path.baseName
            }
            catch( NoSuchFileException e ) {
                throw new IllegalArgumentException("Spack environment file does not exist: $spackEnv")
            }
            catch( Exception e ) {
                throw new IllegalArgumentException("Error parsing Spack environment YAML file: $spackEnv -- Check the log file for details", e)
            }
        }
        // it's interpreted as user provided prefix directory
        else if( spackEnv.contains('/') ) {
            final prefix = spackEnv as Path
            if( !prefix.isDirectory() )
                throw new IllegalArgumentException("Spack prefix path does not exist or is not a directory: $prefix")
            if( prefix.fileSystem != FileSystems.default )
                throw new IllegalArgumentException("Spack prefix path must be a POSIX file path: $prefix")

            return prefix
        }
        else if( spackEnv.contains('\n') ) {
            throw new IllegalArgumentException("Invalid Spack environment definition: $spackEnv")
        }
        else {
            content = spackEnv
        }

        if( arch ) content += arch

        final hash = CacheHelper.hasher(content).hash().toString()
        getCacheDir().resolve("$name-$hash")
    }

    /**
     * Run the spack tool to create an environment in the file system.
     *
     * @param spackEnv The spack environment definition
     * @return the spack environment prefix {@link Path}
     */
    @PackageScope
    Path createLocalSpackEnv(String spackEnv, String arch) {
        final prefixPath = spackPrefixPath(spackEnv, arch)

        final file = new File("${prefixPath.parent}/.${prefixPath.name}.lock")
        final wait = "Another Nextflow instance is creating the spack environment $spackEnv -- please wait till it completes"
        final err =  "Unable to acquire exclusive lock after $createTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: createTimeout, waitMessage: wait, errorMessage: err)

        if( prefixPath.isDirectory() ) {
            log.debug "spack found local env for environment=$spackEnv; path=$prefixPath"
            try {
                // have to check environment packages, because they are stored externally in the host
                // and might get modified there
                mutex .lock { checkLocalSpackEnv0(spackEnv, prefixPath) }
            }
            finally {
                file.delete()
            }
            return prefixPath
        }

        try {
            mutex .lock { createLocalSpackEnv0(spackEnv, prefixPath, arch) }
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
    Path checkLocalSpackEnv0(String spackEnv, Path prefixPath) {

        log.debug "Checking env using spack: $spackEnv [cache $prefixPath]"

        String opts = checksum ? '' : "-n "
        opts += parallelBuilds ? "-j $parallelBuilds " : ''
        opts += '-y '

        def cmd
        cmd =  "spack env activate ${Escape.path(prefixPath)} ; "
        cmd += "spack install ${opts}; "
        cmd += "spack env deactivate"

        try {
            runCommand( cmd )
            log.debug "'spack' check complete env=$spackEnv path=$prefixPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted image file
            prefixPath.delete()
            throw e
        }
        return prefixPath
    }

    @PackageScope
    Path createLocalSpackEnv0(String spackEnv, Path prefixPath, String arch) {

        log.info "Creating env using spack: $spackEnv [cache $prefixPath]"

        String opts = checksum ? '' : "-n "
        opts += parallelBuilds ? "-j $parallelBuilds " : ''
        opts += '-y '

        def cmd
        if( isYamlFilePath(spackEnv) ) {
            cmd =  "spack env create -d ${Escape.path(prefixPath)} ${Escape.path(makeAbsolute(spackEnv))} ; "
            cmd += "spack env activate ${Escape.path(prefixPath)} ; "
            cmd += "spack env view enable ; "
            cmd += "spack config add concretizer:unify:true ; "
            cmd += "spack config add concretizer:reuse:false ; "
            if( arch ) cmd += "spack config add packages:all:target:[$arch] ; "
            cmd += "spack concretize -f ; "
            cmd += "spack install ${opts}; "
            cmd += "spack env deactivate"
        }

        else {
            cmd =  "spack env create -d ${Escape.path(prefixPath)} ; "
            cmd += "spack env activate ${Escape.path(prefixPath)} ; "
            cmd += "spack add $spackEnv ; "
            cmd += "spack env view enable ; "
            cmd += "spack config add concretizer:unify:true ; "
            cmd += "spack config add concretizer:reuse:false ; "
            if( arch ) cmd += "spack config add packages:all:target:[$arch] ; "
            cmd += "spack concretize -f ; "
            cmd += "spack install ${opts}; "
            cmd += "spack env deactivate"
        }

        try {
            runCommand( cmd )
            log.debug "'spack' create complete env=$spackEnv path=$prefixPath"
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
        log.trace """spack env create
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
            def msg = "Failed to create spack environment\n  command: $cmd\n  status : $status\n  message:\n"
            msg += err.toString().trim().indent('    ')
            throw new IllegalStateException(msg)
        }
    }

    /**
     * Given a remote image URL returns a {@link DataflowVariable} which holds
     * the local image path.
     *
     * This method synchronise multiple concurrent requests so that only one
     * image download is actually executed.
     *
     * @param spackEnv
     *      Spack environment string
     * @return
     *      The {@link DataflowVariable} which hold (and pull) the local image file
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String spackEnv, String arch) {

        def spackEnvArch

        spackEnvArch = arch ? "${spackEnv}_${arch}".toString() : spackEnv

        if( spackEnv in spackPrefixPaths ) {
            def msg = "spack found local environment `$spackEnv`"
            if( arch ) msg += " and arch `$arch`"
            log.trace "$msg"
            return spackPrefixPaths[spackEnvArch]
        }

        synchronized (spackPrefixPaths) {
            def result = spackPrefixPaths[spackEnvArch]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ createLocalSpackEnv(spackEnv, arch) })
                spackPrefixPaths[spackEnvArch] = result
            }
            else {
                def msg = "spack found local cache for environment `$spackEnv` (2)"
                if( arch ) msg += " and arch `$arch`"
                log.trace "$msg"
            }
            return result
        }
    }

    /**
     * Create a spack environment caching it in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * environment is actually created.
     *
     * @param spackEnv The spack environment string
     * @return the local environment path prefix {@link Path}
     */
    Path getCachePathFor(String spackEnv, String arch) {
        def promise = getLazyImagePath(spackEnv, arch)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot create Spack environment `$spackEnv`")
        log.trace "Spack cache for env `$spackEnv` path=$result"
        return result
    }

}
