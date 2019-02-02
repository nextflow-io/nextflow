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

package nextflow.file

import java.nio.channels.FileLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.Duration

/**
 * Implements a mutual exclusive file lock strategy to prevent
 * two or more JVMs/process to access the same file or resource
 *
 * NOTE: it should NOT be used to synchronise concurrent from different
 * threads in the same JVM
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileMutex {

    String waitMessage

    String errorMessage

    Duration timeout

    File target

    long delayMillis = 100

    private long timeoutMillis

    FileMutex() { }

    FileMutex(File target) {
        this.target = target
    }

    FileMutex setTimeout( value ) {
        if( value instanceof Duration )
            timeoutMillis = value.millis
        else if( value instanceof CharSequence )
            timeoutMillis = Duration.of(value.toString()).millis
        else if( value instanceof Number )
            timeoutMillis = value.longValue()
        else if( value != null )
            throw new IllegalArgumentException("Not a valid Duration value: $value [${value.class.name}]")

        return this
    }


    FileMutex setWaitMessage( String str ) {
        this.waitMessage = str
        return this
    }

    FileMutex setErrorMessage( String str ) {
        this.errorMessage = str
        return this
    }


    def <T> T lock(Closure<T> closure) {
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
                if( waitMessage && !showed ) {
                    log.info waitMessage
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
