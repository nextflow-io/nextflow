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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.FileVisitResult
import java.nio.file.FileVisitor
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import nextflow.data.cid.serde.CidEncoder
import nextflow.data.cid.serde.CidSerializable
import nextflow.data.config.DataConfig
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.util.TestOnly

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
    private CidEncoder encoder


    DefaultCidStore open(DataConfig config) {
        location = toLocationPath(config.store.location)
        metaLocation = location.resolve(METADATA_PATH)
        encoder = new CidEncoder()
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
            return encoder.decode(path.text) as CidSerializable
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
    void close() throws IOException { }

    @Override
    List<CidSerializable> search(String queryString) {

        def params = null
        if (queryString) {
            params = CidUtils.parseQuery(queryString)
        }
        return searchAllFiles(params)

    }

    private List<CidSerializable> searchAllFiles (Map<String,String> params) {
        final results = new LinkedList<CidSerializable>()

        Files.walkFileTree(metaLocation, new FileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                if (file.name.startsWith('.data.json') ) {
                    final cidObject = encoder.decode(file.text)
                    if (CidUtils.checkParams(cidObject, params)){
                        results.add(cidObject as CidSerializable)
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
