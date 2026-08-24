/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.serde.LinSerializable
/**
 * Default Implementation for the a lineage store.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DefaultLinStore implements LinStore {

    private static String HISTORY_FILE_NAME = ".history"
    private static final String METADATA_FILE = '.data.json'
    private static final String DEFAULT_LOCATION = '.lineage'

    private Path location
    private LinHistoryLog historyLog
    private LinEncoder encoder

    DefaultLinStore open(LineageConfig config) {
        location = toLocationPath(config.store.location)
        encoder = new LinEncoder()
        if( !Files.exists(location) && !Files.createDirectories(location) ) {
            throw new AbortOperationException("Unable to create lineage store directory: $location")
        }
        historyLog = new DefaultLinHistoryLog(location.resolve(HISTORY_FILE_NAME))
        return this
    }

    protected Path toLocationPath(String location) {
        return location
            ? FileHelper.toCanonicalPath(location)
            : Path.of('.').toAbsolutePath().normalize().resolve(DEFAULT_LOCATION)
    }

    @Override
    void save(String key, LinSerializable value) {
        final path = location.resolve("$key/$METADATA_FILE")
        Files.createDirectories(path.parent)
        log.debug "Save LID file path: $path"
        path.text = encoder.encode(value)
    }

    @Override
    LinSerializable load(String key) {
        final path = location.resolve("$key/$METADATA_FILE")
        log.debug("Loading from path $path")
        if( path.exists() )
            return encoder.decode(path.text) as LinSerializable
        log.debug("File for key $key not found")
        return null
    }

    Path getLocation() {
        return location
    }

    @Override
    LinHistoryLog getHistoryLog() {
        return historyLog
    }

    @Override
    void close() throws IOException {}

    @Override
    Stream<String> search(Map<String, List<String>> params) {
        return Files.walk(location)
            .filter { Path path ->
                Files.isRegularFile(path) && path.fileName.toString().startsWith('.data.json')
            }
            .map { Path path ->
                final obj = encoder.decode(path.text)
                final key = location.relativize(path.parent).toString()
                return new AbstractMap.SimpleEntry<String, LinSerializable>(key, obj)
            }
            .filter { entry ->
                LinUtils.checkParams(entry.value, params)
            }
            .map {it->  it.key }
    }

    @Override
    Stream<String> getSubKeys(String parentKey) {
        final startPath = location.resolve(parentKey)

        return Files.walk(startPath)
            .filter { Path path ->
                Files.isRegularFile(path) && path.fileName.toString().startsWith('.data.json') && path.parent != startPath
            }
            .map { Path path ->
                location.relativize(path.parent).toString()
            }
    }
}
