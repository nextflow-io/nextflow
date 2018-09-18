/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.splitter

import java.nio.file.Path

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
        KryoHelper.serialize(allPaths, marker)
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
