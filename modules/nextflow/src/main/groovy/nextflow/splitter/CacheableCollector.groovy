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

package nextflow.splitter

import java.nio.file.AtomicMoveNotSupportedException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.util.KryoHelper
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Defines a trait which can cache the produced chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait CacheableCollector  {

    static final private Logger log = LoggerFactory.getLogger(CacheableCollector)

    @PackageScope Path baseFile

    @PackageScope HashCode hashCode

    @PackageScope List<Path> allPaths = []

    void markComplete() throws IOException {
        // save the list of all chunks
        def marker = getMarkerFile()
        log.trace "Caching chunk paths > marker=$marker; chunks=$allPaths"
        // The marker may be shared by concurrent splitter calls on the same input, and
        // `checkCached` reads it without any locking. Write it to a temporary file and
        // move it into place, so that a reader never sees a partially written marker.
        def temp = marker.resolveSibling("${marker.name}.${UUID.randomUUID()}.tmp".toString())
        try {
            KryoHelper.serialize(allPaths, temp)
            try {
                Files.move(temp, marker, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING)
            }
            catch( AtomicMoveNotSupportedException e ) {
                // file system without atomic rename (e.g. some object stores)
                Files.move(temp, marker, StandardCopyOption.REPLACE_EXISTING)
            }
        }
        finally {
            Files.deleteIfExists(temp)
        }
    }

    Path getMarkerFile() {
        def fileName = ".chunks.${baseFile.name}"
        if( hashCode ) fileName += ".$hashCode"
        return baseFile.resolveSibling(fileName)
    }

    Path getBaseFile() { baseFile }

    boolean checkCached() {
        try {
            def marker = getMarkerFile()
            allPaths = (List<Path>)KryoHelper.deserialize(marker)
            log.trace "Found cached chunk paths > marker=$marker; chunks=$allPaths"
            return true
        }
        catch( IOException e ) {
            return false
        }
    }

    List<Path> getAllChunks() {
        new ArrayList<Path>(allPaths)
    }
}
