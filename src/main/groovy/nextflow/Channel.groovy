package nextflow

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.ControlMessage
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.util.Duration
import nextflow.util.FileHelper
/**
 * Channel factory object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Channel {

    static ControlMessage STOP = PoisonPill.getInstance()

    /**
     * Create an empty channel
     *
     * @return
     */
    static <T> DataflowQueue<T> create() { new DataflowQueue() }

    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */
    static <T> DataflowQueue<T> from( Collection<T> items) { Nextflow.channel(items) }


    /**
     * Creates a channel sending the items in the collection over it
     *
     * @param items
     * @return
     */

    static <T> DataflowQueue<T> from( T... items) { Nextflow.<T>channel(items) }

    /**
     * Convert an object into a *channel* variable that emits that object
     *
     * @param obj
     * @return
     */
    static <T> DataflowVariable just( T obj ) {

        def result = new DataflowVariable<T>()
        result.bind(obj)
        return result

    }

    /**
     * Creates a channel emitting a sequence of integers spaced by a given time interval
     *
     * @param duration
     * @return
     */
    static DataflowQueue interval(String duration) {

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

    static DataflowQueue interval(String duration, Closure closure ) {

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

    static DataflowQueue<Path> files( String filePattern ) {
        assert filePattern

        def channel = create()
        int p = filePattern.indexOf('*')

        if( p == -1 ) {
            channel << FileHelper.asPath( filePattern )
            return channel
        }


        // split the folder and the pattern
        def ( String folder, String pattern ) = getFolderAndPattern(filePattern)

        // now apply glob file search
        def path = FileHelper.asPath(folder)
        def glob = "glob:${folder}${pattern}"
        def matcher = path.getFileSystem().getPathMatcher(glob)

        Thread.start {

            Files.walkFileTree(path, new SimpleFileVisitor<Path>() {
                @Override
                public FileVisitResult visitFile(Path file, BasicFileAttributes attr) throws IOException {
                    if (matcher.matches(file)) { channel.bind(file) }
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

        def folder = null
        def pattern = null
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

}