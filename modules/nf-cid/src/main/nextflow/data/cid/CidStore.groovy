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

import groovy.transform.CompileStatic
import nextflow.data.cid.serde.CidSerializable
import nextflow.data.config.DataConfig
/**
 * Interface for the CID store
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface CidStore extends Closeable {

    /**
     * Open the CID store.
     * @param config Configuration to open the CID store.
     */
    CidStore open(DataConfig config)

    /**
     * Save a CID entry in the store for in a given key.
     * @param key Entry key.
     * @param value Entry object.
     */
    void save(String key, CidSerializable value)

    /**
     * Load an entry for a given CID key.
     * @param key CID key.
     * @return entry value, or null if key does not exists
     */
     CidSerializable load(String key)

    /**
     * Get the {@link CidHistoryLog} object associated to the CidStore.
     * @return {@link CidHistoryLog} object
     */
    CidHistoryLog getHistoryLog()

    List<Object> query(URI queryString)

}
