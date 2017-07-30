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
import nextflow.exception.ProcessNotRecoverableException
import nextflow.util.Escape
/**
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

    SingularityCache(ContainerConfig config, Map<String,String> env=null) {
        this.config = config
        this.env = env ?: System.getenv()
    }

    @PackageScope
    String simpleName(String imageUrl) {
        def p = imageUrl.indexOf('://')
        def name = p != -1 ? imageUrl.substring(p+3) : imageUrl
        name = name.replace(':','-').replace('/','-')
        name += ".img"
        return name
    }

    @PackageScope
    Path checkPath(String str) {
        def result = Paths.get(str)
        if( !result.exists() && !result.mkdirs() ) {
            throw new IOException("Failed to create Singularity cache directory: $str -- Make sure a file with the same does not exist and you have write permission")
        }
        return result
    }

    @PackageScope
    Path getCacheDir() {

        def str = config.cacheDir as String
        if( str )
            return checkPath(str)

        str = env.get('NXF_SINGULARITY_CACHEDIR')
        if( str )
            return checkPath(str)

        def workDir = Global.session.workDir
        if( workDir.fileSystem != FileSystems.default ) {
            throw new IOException("Cannot store Singularity image to a remote work directory -- Use a POSIX compatible work directory or specify an alternative path with the `NXF_SINGULARITY_CACHEDIR` env variable")
        }

        missingCacheDir = true
        def result = workDir.resolve('singularity')
        result.mkdirs()

        return result
    }

    @PackageScope
    Path localImageName(String imageUrl) {
        getCacheDir().resolve( simpleName(imageUrl) )
    }

    @PackageScope
    Path downloadSingularityImage(String imageUrl) {
        def localPath = localImageName(imageUrl)
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

        String cmd = "singularity pull --name ${Escape.path(localPath)} $imageUrl > /dev/null"
        try {
            runCommand( cmd )
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
    int runCommand( String cmd ) {
        log.trace "Singularity pull command: $cmd"

        def proc = new ProcessBuilder(['bash','-c',cmd]).start()
        proc.waitForOrKill(10 * 60 * 1_000)
        def status = proc.exitValue()
        if( status != 0 ) {
            def msg = "Failed to pull singularity image\n  command: $cmd\n  message:\n"
            msg += proc.err?.text?.indent('  ')
            throw new IllegalStateException(msg)
        }
        return status
    }

    @PackageScope
    DataflowVariable<Path> getLocalImageFileName(String imageUrl) {
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


    Path getCachePathFor(String url) {
        def promise = getLocalImageFileName(url)
        def result = promise.getVal()
        if( !result )
            throw new ProcessNotRecoverableException("Missing singularity file for image `$url`")
        if( promise.isError() )
            throw new ProcessNotRecoverableException(promise.getError())
        log.trace "Singularity cache for `$url` path=$result"
        return result
    }

}
