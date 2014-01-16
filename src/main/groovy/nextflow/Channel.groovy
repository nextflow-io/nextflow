package nextflow
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.util.regex.Pattern

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowChannel
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.util.Duration
import nextflow.util.FileHelper
import org.codehaus.groovy.runtime.NullObject

/**
 * Channel factory object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Channel {

    static ControlMessage STOP = PoisonPill.getInstance()

    static NullObject VOID = NullObject.getNullObject()

    /**
     * Create an empty channel
     *
     * @return
     */
    static <T> DataflowChannel<T> create() { new DataflowQueue() }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    static <T> DataflowChannel<T> from( Collection<T> items ) { Nextflow.channel(items) }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */

    static <T> DataflowChannel<T> from( T... items ) { Nextflow.<T>channel(items) }

    /**
     * Convert an object into a *channel* variable that emits that object
     *
     * @param obj
     * @return
     */
    static <T> DataflowChannel<T> just( T obj = null ) {

        def result = new DataflowVariable<T>()
        if( obj != null ) result.bind(obj)

        return result

    }

    /**
     * Creates a channel emitting a sequence of integers spaced by a given time interval
     *
     * @param duration
     * @return
     */
    static DataflowChannel interval(String duration) {

        interval( duration, { index -> index })

    }

    /**
     * Creates a channel emitting a sequence of value given by the closure and spaced by a given time interval.
     *
     * To stop the interval return the special value {@code #STOP}
     *
     * @param duration
     * @return
     */

    static DataflowChannel interval(String duration, Closure closure ) {

        def millis = Duration.create(duration).toMillis()
        def timer = new Timer()
        def channel = create()
        long index = 0

        def task = {
            def value = closure.call(index++)
            channel << value
            if( value == STOP ) {
                timer.cancel()
            }
        }

        timer.schedule( task as TimerTask, millis )

        return channel
    }

    static DataflowChannel<Path> files( Pattern filePattern ) {
        assert filePattern
        // split the folder and the pattern
        def ( String folder, String pattern ) = getFolderAndPattern(filePattern.toString())
        filesImpl( 'regex', folder, pattern, false )
    }

    static DataflowChannel<Path> files( String filePattern ) {
        assert filePattern

        boolean glob  = false
        glob |= filePattern.contains('*')
        glob |= filePattern.contains('?')
        glob = glob ||  filePattern ==~ /.*\{.+\,.+\}.*/

        if( !glob ) {
            return just( FileHelper.asPath(filePattern) )
        }


        // split the folder and the pattern
        def ( String folder, String pattern ) = getFolderAndPattern(filePattern)

        filesImpl('glob', folder, pattern, pattern.startsWith('*')  )
    }

    static private DataflowChannel<Path> filesImpl( String syntax, String folder, String pattern, boolean skipHidden )  {
        assert syntax in ['regex','glob']
        log.debug "files for syntax: $syntax; folder: $folder; pattern: $pattern; skipHidden: $skipHidden"

        // now apply glob file search
        def path = FileHelper.asPath(folder)
        def glob = "$syntax:${folder}${pattern}"
        def matcher = path.getFileSystem().getPathMatcher(glob)
        def channel = new DataflowQueue<Path>()

        Thread.start {

            Files.walkFileTree(path, new SimpleFileVisitor<Path>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attr) throws IOException {
                    if (matcher.matches(file) && ( !skipHidden || !Files.isHidden(file) )) {
                        channel.bind(file)
                    }
                    return FileVisitResult.CONTINUE;
                }

                @Override
                public FileVisitResult visitFileFailed(Path file, IOException exc) throws IOException {
                    return FileVisitResult.CONTINUE;
                }
            })

            channel << STOP
        }

        return channel
    }

    @PackageScope
    static List<String> getFolderAndPattern( String filePattern ) {

        def scheme = null;
        int i = filePattern.indexOf('://')
        if( i != -1 ) {
            scheme = filePattern.substring(0, i+3)
            filePattern = filePattern.substring(i+3)
        }

        def folder
        def pattern
        int p = filePattern.indexOf('*')
        if( p != -1 ) {
            i = filePattern.substring(0,p).lastIndexOf('/')
        }
        else {
            i = filePattern.lastIndexOf('/')
        }

        if( i != -1 ) {
            folder = filePattern.substring(0,i+1)
            pattern = filePattern.substring(i+1)
        }
        else {
            folder = '.'
            pattern = filePattern
        }

        if( scheme ) {
            folder = scheme + folder
        }

        return [ folder, pattern ]

    }

    static DataflowChannel readLines( Map options = [:], Object source ) {
        assert source != null
        assert options != null

        options.into = new DataflowQueue()
        return (DataflowQueue) source.chopLines(options)
    }


    static DataflowChannel readLines( Map options = [:], Object source, Closure closure ) {
        assert source != null
        assert options != null

        options.into = new DataflowQueue()
        return (DataflowQueue) source.chopLines(options,closure)
    }

    static DataflowChannel readFasta( Map options = [:], Object source ) {
        assert source != null
        assert options != null

        options.into = new DataflowQueue()
        return (DataflowQueue) source.chopFasta(options)
    }

    static DataflowChannel readFasta( Map options = [:], Object source, Closure closure ) {
        assert source != null
        assert options != null

        options.into = new DataflowQueue()
        return (DataflowQueue) source.chopFasta(options, closure)
    }


}