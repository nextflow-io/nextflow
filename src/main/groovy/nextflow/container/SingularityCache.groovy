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

    /** Only for debugging purpose - do not use */
    @PackageScope
    SingularityCache() {}

    /**
     * Create a Singularity cache object
     *
     * @param config A {@link ContainerConfig} object
     * @param env The environment configuration object. Specifying {@code null} the current system environment is used
     */
    SingularityCache(ContainerConfig config, Map<String,String> env=null) {
        this.config = config
        this.env = env ?: System.getenv()
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
        name = name.replace(':','-').replace('/','-')
        name += ".img"
        return name
    }

    /**
     * Create the specified directory if not exists
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
            throw new IOException("Failed to create Singularity cache directory: $str -- Make sure a file with the same does not exist and you have write permission")
        }
        return result.toAbsolutePath()
    }

    /**
     * Retried the directory where store the singularity images once downloaded.
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

        def str = config.cacheDir as String
        if( str )
            return checkDir(str)

        str = env.get('NXF_SINGULARITY_CACHEDIR')
        if( str )
            return checkDir(str)

        str = env.get('SINGULARITY_PULLFOLDER')
        if( str )
            return checkDir(str)

        def workDir = Global.session.workDir
        if( workDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Singularity image to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `NXF_SINGULARITY_CACHEDIR` env variable")
        }

        missingCacheDir = true
        def result = workDir.resolve('singularity')
        result.mkdirs()

        return result
    }

    /**
     * Get the path on the file system where store a remote singularity image
     *
     * @param imageUrl The singularity remote URL
     * @return the container image local {@link Path}
     */
    @PackageScope
    Path localImagePath(String imageUrl) {
        getCacheDir().resolve( simpleName(imageUrl) )
    }

    /**
     * Run the singularity tool to pull a remote image and store in the file system.
     * Requires singularity 2.3.x or later.
     *
     * @param imageUrl The singularity image remote URL
     * @return  the container image local {@link Path}
     */
    @PackageScope
    Path downloadSingularityImage(String imageUrl) {
        def localPath = localImagePath(imageUrl)
        if( localPath.exists() ) {
            log.debug "Singularity found local store for image=$imageUrl; path=$localPath"
            return localPath
        }
        else {
            log.trace "Singularity pulling remote image `$imageUrl`"
        }

        if( missingCacheDir )
            log.warn1 "Singularity cache directory has not been defined -- Remote image will be stored in the path: $localPath.parent"

        log.info "Pulling Singularity image $imageUrl [cache $localPath]"

        String cmd = "singularity pull --name ${Escape.path(localPath.getFileName())} $imageUrl > /dev/null"
        try {
            runCommand( cmd, localPath.parent )
            log.debug "Singularity pull complete image=$imageUrl path=$localPath"
        }
        catch( Exception e ){
            // clean-up to avoid to keep eventually corrupted image file
            localPath.delete()
            throw e
        }
        return localPath
    }

    @PackageScope
    int runCommand( String cmd, Path storePath ) {
        log.trace "Singularity pull command: $cmd"

        def builder = new ProcessBuilder(['bash','-c',cmd])
        // workaround due to Singularity issue --> https://github.com/singularityware/singularity/issues/847#issuecomment-319097420
        builder.environment().put('SINGULARITY_PULLFOLDER', storePath.toString())
        def proc = builder.start()
        proc.waitForOrKill(10 * 60 * 1_000)
        def status = proc.exitValue()
        if( status != 0 ) {
            def msg = "Failed to pull singularity image\n  command: $cmd\n  message:\n"
            msg += proc.err?.text?.trim()?.indent('  ')
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
            log.trace "Singularity found local cache for image `$imageUrl`"
            return localImageNames[imageUrl]
        }

        synchronized (localImageNames) {
            def result = localImageNames[imageUrl]
            if( result == null ) {
                result = new LazyDataflowVariable<Path>({ downloadSingularityImage(imageUrl) })
                localImageNames[imageUrl] = result
            }
            else {
                log.trace "Singularity found local cache for image `$imageUrl` (2)"
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
            throw new IllegalStateException("Cannot pull Singularity image `$url`")
        log.trace "Singularity cache for `$url` path=$result"
        return result
    }

}
