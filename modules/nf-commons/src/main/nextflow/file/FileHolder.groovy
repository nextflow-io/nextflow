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
 */

package nextflow.file


import java.nio.file.FileSystems
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

    final Path origPath

    FileHolder( Path inputFile, Path origPath = null ) {
        assert inputFile
        this.sourceObj = inputFile
        this.storePath = real(inputFile)
        this.stageName = norm(inputFile.getFileName())
        this.origPath = origPath
    }

    FileHolder( def origin, Path path, Path origPath = null ) {
        assert origin != null
        assert path != null

        this.sourceObj = origin
        this.storePath = path
        this.stageName = norm(path.getFileName())
        this.origPath = origPath
    }

    protected FileHolder( def source, Path store, def stageName, Path origPath = null ) {
        this.sourceObj = source
        this.storePath = store
        this.stageName = norm(stageName)
        this.origPath = origPath
    }

    FileHolder withName( def stageName )  {
        new FileHolder( this.sourceObj, this.storePath, stageName )
    }

    Path getStorePath() { storePath }

    String getStageName() { stageName }

    Path getOrigPath() { origPath }

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
            // main reason for this is to resolve symlinks to real file location
            // hence apply only for default file system
            // note: also for Google Cloud storage path it may convert to relative path
            // it may return invalid (relative) paths therefore do not apply it
            return path.getFileSystem() == FileSystems.default ? path.toRealPath() : path
        }
        catch( Exception e ) {
            log.trace "Unable to get real path for: $path"
            return path
        }
    }
}
