package nextflow.file

import java.nio.channels.FileLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.Duration

/**
 * Implements a mutual exclusive file lock strategy to prevent
 * two or more JVMs/process to access the same file or resource
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileLocker {

    String message

    Duration timeout

    File target

    long delayMillis = 100

    private long timeoutMillis

    FileLocker(File target) {
        this.target = target
    }

    FileLocker setTimeout( String str ) {
        this.timeoutMillis = Duration.of(str).millis
        return this
    }

    FileLocker setTimeout( Duration duration ) {
        this.timeoutMillis = duration.millis
        return this
    }

    FileLocker setDuration( long millis ) {
        this.timeoutMillis = millis
        return this
    }

    FileLocker setMessage( String str ) {
        this.message = str
        return this
    }


    def lock(Closure closure) {
        assert target
        final file = new RandomAccessFile(target, "rw")
        try {
            /*
             * wait to acquire a lock
             */
            boolean showed=0
            long max = timeout ? timeout.millis : Long.MAX_VALUE
            long begin = System.currentTimeMillis()
            FileLock lock
            while( !(lock=file.getChannel().tryLock()) )  {
                if( message && !showed ) {
                    log.info message
                    showed = true
                }

                if( System.currentTimeMillis()-begin > max )
                    throw new InterruptedException("Cannot acquire FileLock on $target")

                sleep delayMillis
            }

            /*
             * now it can do the job
             */
            try {
                return closure.call()
            }
            finally {
                lock.close()
            }
        }
        finally {
            file.close()
        }

    }
}
