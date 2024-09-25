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
import nextflow.util.Duration
/**
 * Handle caching of remote Charliecloud images
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Patrick Hüther <patrick.huether@gmail.com>
 * @author Niklas Schandry <niklas@bio.lmu.de>
 */
@Slf4j
@CompileStatic
class CharliecloudCache {

    static final private Map<String,DataflowVariable<Path>> localImageNames = new ConcurrentHashMap<>()

    private ContainerConfig config

    private Map<String,String> env

    private boolean missingCacheDir

    private Duration pullTimeout = Duration.of('20min')

    private String registry

    /** Only for debugging purpose - do not use */
    @PackageScope
    CharliecloudCache() {}

    /**
     * Create a Charliecloud cache object
     *
     * @param config A {@link ContainerConfig} object
     * @param env The environment configuration object. Specifying {@code null} the current system environment is used
     */
    CharliecloudCache(ContainerConfig config, Map<String,String> env=null) {
        this.config = config
        this.env = env ?: System.getenv()
    }

    /**
     * Given a remote image URL string returns a normalised name used to store
     * the image in the local file system
     *
     * @param imageUrl
     *          A Docker remote image url. {@code pditommaso/foo:latest}
     * @return
     *          A file name corresponding to the image URL eg. {@code pditommaso-foo.img}
     */
    @PackageScope
    String simpleName(String imageUrl) {
        def p = imageUrl.indexOf('://')
        def name = p != -1 ? imageUrl.substring(p+3) : imageUrl
        
        // add registry
        if( registry ) {
            if( !registry.endsWith('/') ) {
                registry += '/'
            }
            name = registry + name
        }

        name = name.replace(':','+').replace('/','%')
        return name 
    }

    /**
     * Check if the specified directory exists and if it is valid
     * Charliecloud should take care of creating it (and expects to have done so: https://hpc.github.io/charliecloud/faq.html#storage-directory-seems-invalid)
     *
     * @param
     *      str A path string representing a folder where to store the charliecloud containers once downloaded
     * @return
     *      the directory as a {@link Path} object
     */
    @PackageScope
    Path checkDir(String str) {
        def result = Paths.get(str)
        if( !result.exists() ) {
            log.info "Charliecloud cache directory: $str does not exist -- Charliecloud will attempt to initialize it at the specified location"
        } 
        else if( !result.resolve('img').exists() || !result.resolve('dlcache').exists() ) {
            throw new IOException("Charliecloud cache directory exists but seems invalid: $str -- See https://hpc.github.io/charliecloud/faq.html#storage-directory-seems-invalid")
        }
        return result.toAbsolutePath()
    }

    /**
     * Retrieve the directory where to store the charliecloud images once downloaded.
     * It tries these settings in the following order:
     * 1) If writeFake is enabled, the {@code CH_IMAGE_STORAGE} environment variable.
     * 2) {@code charliecloud.cacheDir} setting in the nextflow config file
     * 3) the {@code NXF_CHARLIECLOUD_CACHEDIR} environment variable
     * 4) the {@code $workDir/charliecloud} path
     *
     * @return
     *      the {@code Path} where store the charliecloud images as flattened directories
     */
    @PackageScope
    Path getCacheDir() {

        if( config.pullTimeout )
            pullTimeout = config.pullTimeout as Duration

        def writeFake = true

        if( config.writeFake ) 
            writeFake = config.writeFake?.toString() == 'true'

        def str = config.cacheDir as String

        def charliecloudImageStorage = env.get('CH_IMAGE_STORAGE')

        if( charliecloudImageStorage && writeFake) {
            return checkDir(charliecloudImageStorage)
        }
            
        if( str ) {
            // If charliecloudImageStorage exists and writeFake is true, we never get here
            if( str.equals( charliecloudImageStorage ) ) {
                throw new Exception("`charliecloud.cacheDir` configuration parameter must be different from env variable `CH_IMAGE_STORAGE`")
            }
            return checkDir(str)
        }

        str = env.get('NXF_CHARLIECLOUD_CACHEDIR')

        if( str ) {
            if( str.equals( charliecloudImageStorage ) ) {
                throw new Exception("`NXF_CHARLIECLOUD_CACHEDIR` env variable must be different from env variable `CH_IMAGE_STORAGE`")
            }
            return checkDir(str)
        }

        def workDir = Global.session.workDir

        if( workDir.fileSystem != FileSystems.default ) {
            throw new IOException("Charliecloud cannot store image in a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `NXF_CHARLIECLOUD_CACHEDIR` env variable")
        }

        missingCacheDir = true
        def result = workDir.resolve('charliecloud')

        return result
    }

    /**
     * Get the path on the file system where a remote charliecloud image is stored at
     *
     * @param imageUrl The charliecloud remote URL
     * @return the container image local {@link Path}
     */
    @PackageScope
    Path localImagePath(String imageUrl) {
        getCacheDir().resolve('img').resolve(simpleName(imageUrl))
    }

    /**
     * Run ch-image to pull a remote image and store it in the file system.
     * Requires charliecloud 0.28 or later.
     *
     * @param imageUrl The docker image remote URL
     * @return  the container image local {@link Path}
     */
    @PackageScope
    Path downloadCharliecloudImage(String imageUrl) {
        final localPath = localImagePath(imageUrl)

        if( localPath.exists() ) {
            log.debug "Charliecloud found local store for image=$imageUrl; path=$localPath"
            return localPath
        }

        int count = 0;
        int maxTries = 5;
        boolean imagePulled = false
        while(!imagePulled) {
            try {
                downloadCharliecloudImage0(imageUrl, localPath)
                imagePulled = true
            } catch (e) {
                if (++count == maxTries) throw e
                log.info "Another image is currently pulled. Attempting again in 30 seconds [$count/$maxTries]"
                Thread.sleep(30000)
            }
        }   
        return localPath

    }


    @PackageScope
    Path downloadCharliecloudImage0(String imageUrl, Path targetPath) {

        if( targetPath.exists() ) {
            // If we're here we're an additional process that has waited for the pulling
            // before we got the mutex to advance here.
            log.debug "Charliecloud found local store for image=$imageUrl; path=$targetPath"
            return targetPath
        }
        log.trace "Charliecloud pulling remote image `$imageUrl`"

        if( missingCacheDir )
            log.warn1 "Charliecloud cache directory has not been defined -- Remote image will be stored in the path: $targetPath.parent.parent -- Use the charliecloud.cacheDir config option or set the NXF_CHARLIECLOUD_CACHEDIR variable to specify a different location"

        
        log.info "Charliecloud pulling image $imageUrl [cache $targetPath]"
            
        String cmd = "ch-image pull -s $targetPath.parent.parent $imageUrl > /dev/null"
        try {
            runCommand( cmd, targetPath )
            log.debug "Charliecloud pull complete image=$imageUrl path=$targetPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted container
            targetPath.deleteDir()
            throw e
        }
        return targetPath
    }

    @PackageScope
    int runCommand( String cmd, Path storePath ) {
        log.trace """Charliecloud pull
                     command: $cmd
                     timeout: $pullTimeout
                     folder : $storePath""".stripIndent(true)

        final max = pullTimeout.toMillis()
        final builder = new ProcessBuilder(['bash','-c',cmd])
        builder.environment().remove('CH_IMAGE_STORAGE')
        final proc = builder.start()
        final err = new StringBuilder()
        final consumer = proc.consumeProcessErrorStream(err)
        proc.waitForOrKill(max)
        def status = proc.exitValue()
        if( status != 0 ) {
            consumer.join()
            def msg = "Charliecloud failed to pull image\n  command: $cmd\n  status : $status\n  hint   : Try and increase charliecloud.pullTimeout in the config (current is \"${pullTimeout}\")\n  message:\n"
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
     * @param imageUrl The Docker image remote URL
     * @return The {@link DataflowVariable} which holds (and pulls) the local image directory
     */
    @PackageScope
    DataflowVariable<Path> getLazyImagePath(String imageUrl) {
        if( imageUrl in localImageNames ) {
            log.trace "Charliecloud found local cache for image `$imageUrl`"
            return localImageNames[imageUrl]
        }

        synchronized (localImageNames) {
            def result = localImageNames[imageUrl]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ downloadCharliecloudImage(imageUrl) })
                localImageNames[imageUrl] = result
            }
            else {
                log.trace "Charliecloud found local cache for image `$imageUrl` (2)"
            }
            return result
        }
    }

    /**
     * Pull remote image, caching the resulting flattened directory in the file system.
     *
     * This method synchronise multiple concurrent requests so that only one
     * image download is actually executed.
     *
     * @param url The docker image remote URL
     * @return the container image local {@link Path}
     */
    Path getCachePathFor(String url) {
        def promise = getLazyImagePath(url)
        def result = promise.getVal()
        if( promise.isError() )
            throw new IllegalStateException(promise.getError())
        if( !result )
            throw new IllegalStateException("Charliecloud cannot pull image `$url`")
        log.trace "Charliecloud cache for `$url` path=$result"
        return result
    }

}