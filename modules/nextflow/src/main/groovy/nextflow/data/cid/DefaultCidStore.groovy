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

import java.nio.file.Files
import java.nio.file.Path
import java.util.function.Consumer

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
    private Path metaLocation
    private Path location

    void open(DataConfig config) {
        location = config.store.location
        metaLocation = getMetadataPath(config)
        if( !Files.exists(metaLocation) && !Files.createDirectories(metaLocation) ) {
            throw new AbortOperationException("Unable to create CID store directory: $metaLocation")
        }
    }

    @Override
    void save(String key, Object value) {
        final path = metaLocation.resolve(key)
        Files.createDirectories(path.parent)
        log.debug "Save CID file path: $path"
        path.text = value
    }

    @Override
    void list(String key, Consumer<String> consumer) {
        for( Path it : Files.walk(metaLocation.resolve(key)) ) {
            final fileKey = metaLocation.relativize(it).toString()
            consumer.accept(fileKey)
        }
    }

    @Override
    Object load(String key) {
        metaLocation.resolve(key).text
    }

    @Override
    Path getPath(){ location }

    @Override
    CidHistoryFile getHistoryFile(){
        return new CidHistoryFile(metaLocation.resolve(HISTORY_FILE_NAME))
    }

    static Path getMetadataPath(DataConfig config){ config.store.location.resolve('.meta') }

}
