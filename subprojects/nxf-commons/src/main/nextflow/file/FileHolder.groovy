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

package nextflow.file
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j

/**
 * Implements a special {@code Path} used to stage files in the work area
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class FileHolder  {

    final def sourceObj

    final Path storePath

    final String stageName

    FileHolder( Path inputFile ) {
        assert inputFile
        this.sourceObj = inputFile
        this.storePath = real(inputFile)
        this.stageName = norm(inputFile.getFileName())
    }

    FileHolder( def origin, Path path ) {
        assert origin != null
        assert path != null

        this.sourceObj = origin
        this.storePath = path
        this.stageName = norm(path.getFileName())
    }

    protected FileHolder( def source, Path store, def stageName ) {
        this.sourceObj = source
        this.storePath = store
        this.stageName = norm(stageName)
    }

    FileHolder withName( def stageName )  {
        new FileHolder( this.sourceObj, this.storePath, stageName )
    }

    @PackageScope
    static FileHolder get( def path, def name = null ) {
        Path storePath = path as Path
        def target = name ? name : storePath.getFileName()
        new FileHolder( path, storePath, target )
    }

    /**
     * Make sure the stage name does not start with a slash character
     *
     * @param path
     * @return The normalised path
     */
    static private String norm(path) {
        def result = path.toString()
        return result.startsWith('/') ? result.substring(1) : result
    }

    static private Path real( Path path ) {
        try {
            return path.toRealPath()
        }
        catch( Exception e ) {
            log.trace "Unable to get real path for: $path"
            return path
        }
    }
}
