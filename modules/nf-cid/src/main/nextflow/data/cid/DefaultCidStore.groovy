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

package nextflow.data.cid

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import nextflow.util.TestOnly

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.data.config.DataConfig
import nextflow.exception.AbortOperationException

/**
 * Default Implementation for the a CID store.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DefaultCidStore implements CidStore {

    private static String HISTORY_FILE_NAME =".history"
    private static final String METADATA_FILE = '.data.json'
    private static final String METADATA_PATH = '.meta'
    private Path metaLocation
    private Path location
    private CidHistoryLog historyLog

    void open(DataConfig config) {
        location = config.store.location
        metaLocation = config.store.location.resolve(METADATA_PATH)
        if( !Files.exists(metaLocation) && !Files.createDirectories(metaLocation) ) {
            throw new AbortOperationException("Unable to create CID store directory: $metaLocation")
        }
        historyLog = new CidHistoryFile(config.store.logLocation ?: metaLocation.resolve(HISTORY_FILE_NAME))
    }

    @Override
    void save(String key, Object value) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        Files.createDirectories(path.parent)
        log.debug "Save CID file path: $path"
        path.text = value
    }

    @Override
    Object load(String key) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        log.debug("Loading from path $path")
        if (path.exists())
            return path.text
        log.debug("File for key $key not found")
        return null
    }

    @Override
    Path getPath(){ location }

    @TestOnly
    Path getMetadataPath() {metaLocation}

    @Override
    CidHistoryLog getHistoryLog(){ historyLog }

}
