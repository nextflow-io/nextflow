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

import groovy.transform.CompileStatic
import nextflow.lineage.serde.LinSerializable
import nextflow.lineage.config.LineageConfig
/**
 * Interface for the lineage store
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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
     * @return Key-lineage entry pairs fulfilling the query params
     */
    Map<String,LinSerializable> search(Map<String, List<String>> params)

}
