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

package nextflow.util

import java.util.concurrent.LinkedBlockingQueue

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * A blocking queue implementation that changes the semantic of
 *  {@link java.util.concurrent.BlockingQueue#offer(java.lang.Object)}
 *  to make behaving in a blocking manner.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BlockingBlockingQueue<E> extends LinkedBlockingQueue<E> {

    BlockingBlockingQueue(int maxSize) {
        super(maxSize);
    }

    /**
     * Overrides the offer method semantic to make it block when
     * the queue reached the max capacity.
     *
     * @param e
     *      The element to be added in the queue
     * @return
     *      {@code true} if the element was added successfully,
     *      {@code false} if the adding operation was interrupted.
     */
    boolean offer(E e) {
        try {
            this.put(e)
            return true
        }
        catch (InterruptedException error) {
            Thread.currentThread().interrupt()
            return false
        }
    }
}