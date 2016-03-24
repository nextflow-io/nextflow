/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import nextflow.util.KryoHelper
/**
 * Defines a trait which can cache the produced chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait CacheableCollector  {

    @PackageScope Path baseFile

    @PackageScope List<Path> allPaths = []

    void markComplete() throws IOException {
        // save the list of all chunks
        def marker = baseFile.resolveSibling('.chunks')
        KryoHelper.serialize(allPaths, marker)
    }

    boolean checkCached() {
        try {
            def marker = baseFile.resolveSibling('.chunks')
            allPaths = (List<Path>)KryoHelper.deserialize(marker)
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
