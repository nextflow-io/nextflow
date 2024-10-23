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

package nextflow.container

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.LazyDataflowVariable
import nextflow.Const
import nextflow.Global
import nextflow.SysEnv
import nextflow.file.FileMutex
import nextflow.util.Duration
import nextflow.util.Escape
/**
 * Handle caching of remote Singularity images
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SingularityCache {

    static final private Map<String,DataflowVariable<Path>> localImageNames = new ConcurrentHashMap<>()

    private ContainerConfig config

    private Map<String,String> env

    private boolean missingCacheDir

    private Duration pullTimeout = Duration.of('20min')

    /** Only for debugging purpose - do not use */
    @PackageScope
    SingularityCache() {}

    protected String getBinaryName() { return 'singularity' }

    protected String getAppName() { getBinaryName().capitalize() }

    protected String getEnvPrefix() { getBinaryName().toUpperCase() }

    /**
     * Create a Singularity cache object
     *
     * @param config A {@link ContainerConfig} object
     * @param env The environment configuration object. Specifying {@code null} the current system environment is used
     */
    SingularityCache(ContainerConfig config, Map<String,String> env=null) {
        this.config = config
        this.env = env ?: SysEnv.get()
    }

    /**
     * Given a remote image URL string returns a normalised name used to store
     * the image in the local file system
     *
     * @param imageUrl
     *          A Singularity remote image url. It must be prefixed with a pseudo-protocol supported by
     *          by the underlying singularity tool eg. {@code docker://pditommaso/foo:latest}
     * @return
     *          A file name corresponding to the image URL eg. {@code pditommaso-foo.img}
     */
    @PackageScope
    String simpleName(String imageUrl) {
        def p = imageUrl.indexOf('://')
        def name = p != -1 ? imageUrl.substring(p+3) : imageUrl
        String extension = '.img'
        if( name.contains('.sif:') ) {
            extension = '.sif'
            name = name.replace('.sif:','-')
        }
        else if( name.endsWith('.sif') ) {
            extension = '.sif'
            name = name.substring(0,name.length()-4)
        }
        name = name.replace(':','-').replace('/','-')
        return name + extension
    }

    /**
     * Create the specified directory if it does not exist
     *
     * @param
     *      str A path string representing a folder where store the singularity images once downloaded
     * @return
     *      the directory as a {@link Path} object
     */
    @PackageScope
    Path checkDir(String str) {
        def result = Paths.get(str)
        if( !result.exists() && !result.mkdirs() ) {
            throw new IOException("Failed to create ${appName} cache directory: $str -- Make sure a file with the same name does not exist and you have write permission")
        }
        return result.toAbsolutePath()
    }

    @PackageScope
    Path existsDir(String str) {
        def result = Paths.get(str)
        if( !result.exists() ) {
            throw new IOException("Missing ${appName} library directory: $str")
        }
        return result.toAbsolutePath()
    }

    /**
     * Retrieve the directory where store the singularity images once downloaded.
     * If tries these setting in the following order:
     * 1) {@code singularity.cacheDir} setting in the nextflow config file;
     * 2) the {@code NXF_SINGULARITY_CACHEDIR} environment variable
     * 3) the {@code $workDir/singularity} path
     *
     * @return
     *      the {@code Path} where store the singularity images
     */
    @PackageScope
    Path getCacheDir() {

        if( config.pullTimeout )
            pullTimeout = config.pullTimeout as Duration

        def str = config.cacheDir as String
        if( str )
            return checkDir(str)

        str = env.get("NXF_${envPrefix}_CACHEDIR".toString())
        if( str )
            return checkDir(str)

        str = env.get("${envPrefix}_PULLFOLDER".toString())
        if( str )
            return checkDir(str)

        def workDir = Global.session.workDir
        if( workDir.fileSystem != FileSystems.default ) {
            // when the work dir is a remote path use the local launch directory to cache image files
            workDir = Const.appCacheDir.toAbsolutePath()
        }

        missingCacheDir = true
        def result = workDir.resolve('singularity')
        result.mkdirs()

        return result
    }

    /**
     * Defines the Singularity *library* path. The library directory is checked by Nextflow
     * before the cache directory to retrieve a image files. This can be useful to provide
     * a read-only catalog of images, and still have the ability to download and cache missing
     * images into a separate (caching) path.
     *
     * The library directory defined using the setting below in the following order:
     * 1) {@code singularity.libraryDir} setting in the nextflow config file;
     * 2) the {@code NXF_SINGULARITY_LIBRARYDIR} environment variable
     *
     * @return The library directory or {@code null} if not defined
     */
    @PackageScope
    Path getLibraryDir() {
        def str = config.libraryDir as String ?: env.get("NXF_${envPrefix}_LIBRARYDIR".toString())
        if( str )
            return existsDir(str)

        return null
    }

    @PackageScope
    Path localLibraryPath(String imageUrl) {
        getLibraryDir()?.resolve( simpleName(imageUrl) )
    }

    /**
     * Get the path on the file system where store a remote singularity image
     *
     * @param imageUrl The singularity remote URL
     * @return the container image local {@link Path}
     */
    @PackageScope
    Path localCachePath(String imageUrl) {
        getCacheDir().resolve( simpleName(imageUrl) )
    }

    @PackageScope
    Path getTempImagePath(Path targetPath) {
        targetPath.resolveSibling("${targetPath.name}.pulling.${System.currentTimeMillis()}")
    }

    /**
     * Run the singularity tool to pull a remote image and store in the file system.
     * Requires singularity 2.3.x or later.
     *
     * @param imageUrl The singularity image remote URL
     * @return  the container image local {@link Path}
     */
    @PackageScope
    Path downloadContainerImage(String imageUrl) {
        // check for the image in the local library dir
        // see https://github.com/nextflow-io/nextflow/issues/1879
        final libraryPath = localLibraryPath(imageUrl)
        if( libraryPath?.exists() ) {
            log.debug "${appName} found local library for image=$imageUrl; path=$libraryPath"
            return libraryPath
        }

        // check for the image in the cache dir
        // if the image does not exist in the cache dir, download it
        final localPath = localCachePath(imageUrl)
        if( localPath.exists() ) {
            log.debug "${appName} found local store for image=$imageUrl; path=$localPath"
            return localPath
        }

        final file = new File("${localPath.parent}/.${localPath.name}.lock")
        final wait = "Another Nextflow instance is pulling the ${appName} image $imageUrl -- please wait the download completes"
        final err =  "Unable to acquire exclusive lock after $pullTimeout on file: $file"

        final mutex = new FileMutex(target: file, timeout: pullTimeout, waitMessage: wait, errorMessage: err)
        try {
            mutex .lock { downloadContainerImage0(imageUrl, localPath) }
        }
        finally {
            file.delete()
        }

        return localPath
    }


    @PackageScope
    Path downloadContainerImage0(String imageUrl, Path targetPath) {

        if( targetPath.exists() ) {
            // If we're here we're an additional process that has waited for the pulling
            // before we got the mutex to advance here.
            log.debug "${appName} found local store for image=$imageUrl; path=$targetPath"
            return targetPath
        }
        log.trace "${appName} pulling remote image `$imageUrl`"

        if( missingCacheDir )
            log.warn1 "${appName} cache directory has not been defined -- Remote image will be stored in the path: $targetPath.parent -- Use the environment variable NXF_${envPrefix}_CACHEDIR to specify a different location"

        log.info "Pulling ${appName} image $imageUrl [cache $targetPath]"

        // Construct a temporary name for the image file
        final tmpFile = getTempImagePath(targetPath)
        final noHttpsOption = (config.noHttps)? '--nohttps' : ''

        String cmd = "${binaryName} pull ${noHttpsOption} --name ${Escape.path(tmpFile.name)} $imageUrl > /dev/null"
        try {
            runCommand( cmd, tmpFile.parent )
            Files.move( tmpFile, targetPath )
            log.debug "${appName} pull complete image=$imageUrl path=$targetPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted image file
            tmpFile.delete()
            throw e
        }
        return targetPath
    }

    @PackageScope
    int runCommand( String cmd, Path storePath ) {
        log.trace """${appName} pull
                     command: $cmd
                     timeout: $pullTimeout
                     folder : $storePath""".stripIndent(true)

        final max = pullTimeout.toMillis()
        final builder = new ProcessBuilder(['bash','-c',cmd])
        // workaround due to Singularity issue --> https://github.com/singularityware/singularity/issues/847#issuecomment-319097420
        builder.directory(storePath.toFile())
        builder.environment().remove("${envPrefix}_PULLFOLDER".toString())
        final proc = builder.start()
        final err = new StringBuilder()
        final consumer = proc.consumeProcessErrorStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            def msg = "Failed to pull singularity image\n  command: $cmd\n  status : $status\n  hint   : Try and increase ${binaryName}.pullTimeout in the config (current is \"${pullTimeout}\")\n  message:\n"
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
     * @param imageUrl The singularity image remote URL
     * @return The {@link DataflowVariable} which hold (and pull) the local image file
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String imageUrl) {
        if( imageUrl in localImageNames ) {
            log.trace "${appName} found local cache for image `$imageUrl`"
            return localImageNames[imageUrl]
        }

        synchronized (localImageNames) {
            def result = localImageNames[imageUrl]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ downloadContainerImage(imageUrl) })
                localImageNames[imageUrl] = result
            }
            else {
                log.trace "${appName} found local cache for image `$imageUrl` (2)"
            }
            return result
        }
    }

    /**
     * Pull a Singularity remote image, caching the resulting image in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * image download is actually executed.
     *
     * @param url The singularity image remote URL
     * @return the container image local {@link Path}
     */
    Path getCachePathFor(String url) {
        def promise = getLazyImagePath(url)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Cannot pull ${appName} image `$url`")
        log.trace "${appName} cache for `$url` path=$result"
        return result
    }

}
