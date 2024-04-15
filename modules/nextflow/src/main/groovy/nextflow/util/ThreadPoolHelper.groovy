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
 *
 */

package nextflow.util

import java.util.concurrent.ExecutorService
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import jdk.internal.vm.ThreadContainer

/**
 * Thread pool helpers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class ThreadPoolHelper {

    static void await(ExecutorService pool, Duration maxAwait, String waitMessage, String exitMsg) {
        final max = maxAwait.millis
        final t0 = System.currentTimeMillis()
        // wait for ongoing file transfer to complete
        int count=0
        while( true ) {
            final terminated = pool.awaitTermination(5, TimeUnit.SECONDS)
            if( terminated )
                break

            final delta = System.currentTimeMillis()-t0
            if( delta > max ) {
                log.warn(exitMsg)
                break
            }

            // log to console every 10 minutes (120 * 5 sec)
            if( count % 120 == 0 ) {
                log.info1(format(waitMessage, pool))
            }
            // log to the debug file every minute (12 * 5 sec)
            else if( count % 12 == 0 ) {
                log.debug(format(waitMessage, pool))
            }
            // increment the count
            count++
        }
    }

    static protected int pending(ExecutorService pool) {
        if( pool instanceof ThreadPoolExecutor) {
            final p1 = ((ThreadPoolExecutor)pool)
            return p1.getTaskCount() - p1.getCompletedTaskCount()
        }
        else if( pool instanceof ThreadContainer ) {
            return pool.threadCount()
        }
        return -1
    }

    static protected String format(String message, ExecutorService pool) {
        final pending = pending(pool)
        return pending>=0
                ? String.format(message, pending)
                : message.replaceAll(/\([^)]+\)/, '') // <-- remove text between parentheses
    }
}
