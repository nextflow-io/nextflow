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

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j

import nextflow.file.FileHolder

/**
 * Implements a special {@code Path} used to stage files in the work area
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class RemoteFileHolder extends FileHolder {

    final Path origPath

    RemoteFileHolder( Path inputFile, Path origPath ) {
        super(inputFile)
        this.origPath = origPath
    }

    RemoteFileHolder( def origin, Path path, Path origPath ) {
        super(origin, path)
        this.origPath = origPath
    }

    protected RemoteFileHolder( def source, Path store, def stageName, Path origPath ) {
        super(source, store, stageName)
        this.origPath = origPath
    }

    Path getOrigPath() { origPath }
}
