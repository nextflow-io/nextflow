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

<<<<<<< HEAD
import nextflow.file.FileHelper
=======
import nextflow.data.cid.serde.Encoder
import nextflow.data.cid.serde.JsonEncoder
>>>>>>> bdbe5bdb0 (Add serde interfaces)
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
    private Encoder<String,CidSerializable> encoder

    DefaultCidStore open(DataConfig config) {
        location = toLocationPath(config.store.location)
        metaLocation = location.resolve(METADATA_PATH)
        encoder = new JsonEncoder<CidSerializable>() {}
        if( !Files.exists(metaLocation) && !Files.createDirectories(metaLocation) ) {
            throw new AbortOperationException("Unable to create CID store directory: $metaLocation")
        }
        historyLog = new CidHistoryFile(metaLocation.resolve(HISTORY_FILE_NAME))
        return this
    }

    protected Path toLocationPath(String location) {
        return location
            ? FileHelper.toCanonicalPath(location)
            : Path.of('.').toAbsolutePath().normalize().resolve('data')
    }

    @Override
    void save(String key, CidSerializable value) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        Files.createDirectories(path.parent)
        log.debug "Save CID file path: $path"
        path.text = encoder.encode(value)
    }

    @Override
    CidSerializable load(String key) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        log.debug("Loading from path $path")
        if (path.exists())
            return encoder.decode(path.text)
        log.debug("File for key $key not found")
        return null
    }

    Path getLocation(){
        return location
    }

    @TestOnly
    Path getMetadataPath() {
        return metaLocation
    }

    @Override
    CidHistoryLog getHistoryLog(){
        return historyLog
    }

    @Override
    void close() { }
}
