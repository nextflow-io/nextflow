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
 */

package nextflow.splitter

import groovy.transform.CompileStatic

/**
 * Defines an abstract contract to collect items produced by splitters into chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface CollectorStrategy {

    /**
     * @param item Add an item to the current chunk
     */
    void add( Object item )

    /**
     * @return {@code true} if at least an item as been added
     */
    boolean hasChunk()

    /**
     * @return The chunk made up of the added items
     */
    Object nextChunk()

}
