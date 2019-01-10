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

package nextflow.splitter

import groovy.transform.CompileStatic

/**
 * Count the number of entries by which each chunk is made up
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class EntryCounter {

    private long increment = 1

    /**
     * The number of entries that are expected to made-up a split chunk
     */
    private long size

    /**
     * The current count of parsed entries
     */
    private long current

    /**
     * Whenever this counter is enabled, by default when the size > 1
     * (ie. when a chunk is composed at least by two parsed items)
     */
    private boolean enabled

    /**
     * Create a splitter entries counter
     *
     * @param size The number of entries by which a chunk is expected to be composed
     */
    EntryCounter(long size) {
        assert size>0
        this.size = size
        this.enabled = size > 1
    }

    /**
     * Create a splitter entries counter
     *
     * @param size The number of entries by which a chunk is expected to be composed
     * @param enabled Whenever this counter logic is enabled
     */
    EntryCounter(long size, boolean enabled) {
        assert size>0
        this.size = size
        this.enabled = enabled
    }

    /**
     * Increment the entries count and check if the chunk is complete
     *
     * @return
     */
    boolean isChunckComplete() {
        current += increment
        current >= size
    }

    /**
     * Reset the current count of parsed entries
     */
    void reset() {
        current=0
    }

    /**
     * @param value The value that is used to increment the counter each time a entry is parsed. By default {code 1}.
     */
    void setIncrement(long value) {
        this.increment = value
    }

    /**
     * @return The number of entries that are expected to made-up a split chunk
     */
    long getSize() {
        size
    }

    boolean isEnabled() {
        enabled
    }

    protected getCounter() {
        counter
    }

    protected getIncrement() {
        increment
    }

}
