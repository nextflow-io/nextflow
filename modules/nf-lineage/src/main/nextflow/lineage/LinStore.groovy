/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage

import com.google.common.annotations.Beta
import groovy.transform.CompileStatic
import java.util.stream.Stream
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.serde.LinSerializable
/**
 * Interface for the lineage store
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Beta
@CompileStatic
interface LinStore extends Closeable {

    /**
     * Open the lineage store.
     * @param config Configuration to open the lineage store.
     */
    LinStore open(LineageConfig config)

    /**
     * Save a lineage entry in the store for in a given key.
     * @param key Entry key.
     * @param value Entry object.
     */
    void save(String key, LinSerializable value)

    /**
     * Load an entry for a given Lineage ID key.
     * @param key LID key.
     * @return entry value, or null if key does not exists
     */
     LinSerializable load(String key)

    /**
     * Get the {@link LinHistoryLog} object associated to the lineage store.
     * @return {@link LinHistoryLog} object
     */
    LinHistoryLog getHistoryLog()

    /**
     * Search for lineage entries.
     * @param params Map of query params
     * @return Stream with keys fulfilling the query params
     */
    Stream<String> search(Map<String, List<String>> params)

    /**
     * Search for keys starting with a parent key.
     * For example, if a LinStore contains the following keys: '123abc', '123abc/samples/file1.txt' and '123abc/summary',
     * The execution of the function with parentKey='123abc' will return a stream with '123abc/samples/file1.txt' and '123abc/summary'.
     * Similarly, the execution of the function with parentKey='123abc/samples' will just return '123abc/samples/file1.txt"
     *
     * @param parentKey
     * @return Stream of keys starting with parentKey
     */
    Stream<String> getSubKeys(String parentKey)
}
