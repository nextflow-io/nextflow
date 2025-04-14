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
import groovy.util.logging.Slf4j

import java.nio.file.FileVisitResult
import java.nio.file.FileVisitor
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.serde.LinSerializable
import nextflow.lineage.config.LineageConfig
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.util.TestOnly

/**
 * Default Implementation for the a lineage store.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DefaultLinStore implements LinStore {

    private static String HISTORY_FILE_NAME =".history"
    private static final String METADATA_FILE = '.data.json'
    private static final String METADATA_PATH = '.meta'

    private Path metaLocation
    private Path location
    private LinHistoryLog historyLog
    private LinEncoder encoder


    DefaultLinStore open(LineageConfig config) {
        location = toLocationPath(config.store.location)
        metaLocation = location.resolve(METADATA_PATH)
        encoder = new LinEncoder()
        if( !Files.exists(metaLocation) && !Files.createDirectories(metaLocation) ) {
            throw new AbortOperationException("Unable to create lineage store directory: $metaLocation")
        }
        historyLog = new DefaultLinHistoryLog(metaLocation.resolve(HISTORY_FILE_NAME))
        return this
    }

    protected Path toLocationPath(String location) {
        return location
            ? FileHelper.toCanonicalPath(location)
            : Path.of('.').toAbsolutePath().normalize().resolve('data')
    }

    @Override
    void save(String key, LinSerializable value) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        Files.createDirectories(path.parent)
        log.debug "Save LID file path: $path"
        path.text = encoder.encode(value)
    }

    @Override
    LinSerializable load(String key) {
        final path = metaLocation.resolve("$key/$METADATA_FILE")
        log.debug("Loading from path $path")
        if (path.exists())
            return encoder.decode(path.text) as LinSerializable
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
    LinHistoryLog getHistoryLog(){
        return historyLog
    }

    @Override
    void close() throws IOException { }

    @Override
    Map<String, LinSerializable> search(String queryString) {
        def params = null
        if (queryString) {
            params = LinUtils.parseQuery(queryString)
        }
        return searchAllFiles(params)
    }

    private Map<String, LinSerializable> searchAllFiles (Map<String,String> params) {
        final results = new HashMap<String, LinSerializable>()

        Files.walkFileTree(metaLocation, new FileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                if (file.name.startsWith('.data.json') ) {
                    final lidObject = encoder.decode(file.text)
                    if (LinUtils.checkParams(lidObject, params)){
                        results.put(metaLocation.relativize(file.getParent()).toString(), lidObject as LinSerializable)
                    }
                }
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFileFailed(Path file, IOException exc) throws IOException {
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
                FileVisitResult.CONTINUE
            }
        })

        return results
    }
}
