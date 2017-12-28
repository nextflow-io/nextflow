/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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

