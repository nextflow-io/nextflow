/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.io

import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

import groovy.util.logging.Slf4j

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ByteDumper extends Thread {

    InputStream fInputStream
    boolean fTerminated
    CountDownLatch barrier = new CountDownLatch(1)
    Closure fCallback
    File fInputFile

    ByteDumper(InputStream input0, Closure callback0 ) {
        assert input0

        this.fInputStream = new BufferedInputStream(input0)
        this.fCallback = callback0
        setDaemon(true)
    }

    ByteDumper( File file0, Closure callback0 ) {
        assert file0

        this.fInputFile = file0
        this.fCallback = callback0
        setDaemon(true)
    }

    /**
     *  Interrupt the dumper thread
     */
    def void terminate() { fTerminated = true }

    /**
     * Await that the thread finished to read the process stdout
     *
     * @param millis Maximum time (in millis) to await
     */
    def void await(long millis=0) {
        if( millis ) {
            barrier.await(millis, TimeUnit.MILLISECONDS)
        }
        else {
            barrier.await()

        }
    }


    @Override
    def void run() {

        try {
            consume()
        }
        finally{
            if( fInputStream ) fInputStream.closeQuietly()
            barrier.countDown()
        }

    }

    def void consume() {

        // -- if a file reference has been provided instead of stream
        //    wait for the file to exist
        while( fInputStream == null && fInputFile !=null ) {
            if( fTerminated ) {
                log.trace "consume '${getName()}' -- terminated"
                return
            }

            if( fInputFile.exists() ) {
                log.trace "consume '${getName()}' -- ${fInputFile} exists: true"
                fInputStream = new BufferedInputStream(new FileInputStream(fInputFile))
            }
            else {
                sleep 200
            }
        }

        if ( !fCallback ) {
            log.trace "consume '${getName()}' -- no callback"

            return
        }

        byte[] buf = new byte[8192]
        int next
        try {
            while ((next = fInputStream.read(buf)) != -1 && !fTerminated) {
                log.trace  "consume '${getName()}' -- reading "
                fCallback.call(buf, next)
            }

            log.trace "consume '${getName()}' -- exit -- terminated: ${fTerminated}"

        } catch (IOException e) {
            throw new RuntimeException("exception while dumping process stream", e)
        }
    }
}

