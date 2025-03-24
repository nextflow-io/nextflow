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

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

import java.nio.file.FileVisitResult
import java.nio.file.FileVisitor
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import nextflow.data.cid.fs.CidPath
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
    void close() throws IOException { }

    @Override
    List<Object> query(URI uri) {
        log.debug("Query $uri.query, Path: $uri.path, Scheme: $uri.path $uri.authority $uri.rawPath $uri.rawAuthority")
        def params = null
        if (uri.query) {
            params = uri.query.split('&').collectEntries {
                it.split('=').collect { URLDecoder.decode(it, 'UTF-8') }
            } as Map<String, String>
        }
        String key = uri.authority ? uri.authority + uri.path : uri.path
        if (key == CidPath.SEPARATOR) {
            searchAllFiles(params)
        } else {
            searchPath(key, params)
        }

    }

    private List<Object> searchAllFiles (Map<String,String> params) {
        final results = new LinkedList<Object>()
        final slurper = new JsonSlurper()

        Files.walkFileTree(metaLocation, new FileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                FileVisitResult.CONTINUE
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {

                if (file.name.startsWith('.data.json') ) {
                    final cidObject = slurper.parse( file.text.toCharArray() ) as Map
                    DefaultCidStore.treatObject(cidObject, params, results)
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

    private List<Object> searchPath( String path, Map<String,String> params, String[] childs = []) {
        final results = new LinkedList<Object>()
        final slurper = new JsonSlurper()
        final object = load(path)
        if ( object ) {
            final cidObject = slurper.parse(object.toString().toCharArray()) as Map
            if (childs && childs.size() > 0) {
                final output = cidObject.navigate(childs.join('.'))
                if (output) {
                    treatObject(output, params, results)
                } else {
                    throw new FileNotFoundException("Cid object $path/${childs ? childs.join('/') : ''} not found.")
                }
            } else {
                treatObject(cidObject, params, results)
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(path)
            final parent = currentPath.getParent()
            if( parent) {
                ArrayList<String> newChilds = new ArrayList<String>()
                newChilds.add(currentPath.getFileName().toString())
                newChilds.addAll(childs)
                return searchPath(parent.toString(), params, newChilds as String[])
            } else {
                throw new FileNotFoundException("Cid object $path/${childs ? childs.join('/') :''} not found.")
            }
        }
        return results
    }

    protected static void treatObject(def object, Map<String,String> params, List<Object> results) {
        if (params) {
            if (object instanceof Collection) {
                (object as Collection).forEach { treatObject(it, params, results) }
            } else if (checkParams(object as Map, params)) {
                results.add(object)
            }
        } else {
            results.add(object)
        }
    }

    private static boolean checkParams(Map object, Map<String,String> params) {
        log.debug("Checking $object, $params")
        for (final entry : params.entrySet()) {
            final value = object.navigate(entry.key)
            log.debug("comparing $value, $entry.value")
            if (!value || value.toString() != entry.value.toString() ) {
                return false
            }
        }
        log.debug("Returning true")
        return true
    }

}
