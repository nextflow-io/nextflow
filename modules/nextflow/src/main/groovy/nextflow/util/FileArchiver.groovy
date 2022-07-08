/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.util

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.file.FileHelper

/**
 * Helper class to resolve bucket archive paths
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Singleton(lazy = true, strict = false)
class FileArchiver {

    private final Path baseDir
    private final Path targetDir

    Path getBaseDir() { return baseDir }

    Path getTargetDir() { return targetDir }

    boolean asBoolean() {
        return baseDir!=null && targetDir!=null
    }

    FileArchiver() {
        this(System.getenv())
    }

    FileArchiver(Map<String,String> env) {
        final paths = env.get('NXF_ARCHIVE_DIR')?.tokenize(',')
        if( paths ) {
            if( paths.size()!=2 )
                throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR format - expected exactly two paths separated by a command - offending value: ${System.getenv('NXF_ARCHIVE_DIR')}")
            if( !paths[0].startsWith('/') )
                throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR base path - it must start with a slash character - offending value: '${paths[0]}'")
            final scheme = FileHelper.getUrlProtocol(paths[1])
            if ( !scheme )
                throw new IllegalArgumentException("Invalid NXF_ARCHIVE_DIR target path - it must start be a remote path - offending value: '${paths[1]}'")
            this.baseDir = Path.of(paths[0])
            this.targetDir = FileHelper.asPath(paths[1])
        }
    }

    Path archivePath(Path source) {
        if( baseDir==null )
            return null
        if( source==null )
            return null
        if( !source.startsWith(baseDir) )
            return null
        final delta = baseDir.relativize(source)
        // convert to string to prevent 'ProviderMismatchException'
        return targetDir.resolve(delta.toString())
    }
}
