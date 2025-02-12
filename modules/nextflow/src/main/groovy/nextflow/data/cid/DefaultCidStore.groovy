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
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DefaultCidStore implements CidStore {

    private Path location

    void open(DataConfig config) {
        location = config.store.location.resolve('.meta')
        if( !Files.exists(location) && !Files.createDirectories(location) ) {
            throw new AbortOperationException("Unable to create CID store directory: $location")
        }
    }

    @Override
    void save(String key, Object value) {
        final path = location.resolve(key)
        Files.createDirectories(path.parent)
        log.debug "Save CID file path: $path"
        path.text = value
    }

    @Override
    void list(String key, Consumer<String> consumer) {
        for( Path it : Files.walk(location.resolve(key)) ) {
            final fileKey = location.relativize(it).toString()
            consumer.accept(fileKey)
        }
    }

    @Override
    Object load(String key) {
        location.resolve(key).text
    }

    @Override
    Path getPath(){ location }


}
