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

import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper
import nextflow.util.CacheHelper.HashMode
/**
 * Implements a special {@code Path} used to stage files in the work area
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class FileHolder implements CacheFunnel {

    final def sourceObj

    final Path storePath

    final String stageName

    FileHolder( Path path ) {
        this.sourceObj = path
        this.storePath = real(path)
        this.stageName = norm(path.getFileName())
    }

    FileHolder( def origin, Path path ) {
        this.sourceObj = origin
        this.storePath = real(path)
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

    Path getSourcePath() {
        sourceObj instanceof Path ? sourceObj : null
    }

    Path getStorePath() { storePath }

    String getStageName() { stageName }

    @Override
    Hasher funnel(Hasher hasher, HashMode mode) {
        return CacheHelper.hasher(hasher, sourceObj, mode)
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
