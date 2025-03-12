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
 *
 */

package nextflow.data.cid

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.data.config.DataConfig

/**
 * Interface for the CID store
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface CidStore {

    /**
     * Open the CID store.
     * @param config Configuration to open the CID store.
     */
    void open(DataConfig config)

    /**
     * Save a CID entry in the store for in a given key.
     * @param key Entry key.
     * @param value Entry object.
     */
    void save(String key, Object value)

    /**
     * Load an entry for a given CID key.
     * @param key CID key.
     * @return entry value, or null if key does not exists
     */
    Object load(String key)

    /**
     * Get the CID store location path.
     * @return CID store location path.
     */
    Path getPath()

    /**
     * Get the {@link CidHistoryLog} object associated to the CidStore.
     * @return {@link CidHistoryLog} object
     */
    CidHistoryLog getHistoryLog()

}
