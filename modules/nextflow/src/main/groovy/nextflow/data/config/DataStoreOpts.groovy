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

package nextflow.data.config

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.data.fs.CidFileSystemProvider
import nextflow.file.FileHelper

/**
 * Model data store options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataStoreOpts {

    final Path location

    DataStoreOpts(Map opts) {
        this.location = opts.location
            ? loc(opts.location as String)
            : Path.of('.').toAbsolutePath().normalize().resolve('data')
    }

    static private Path loc(String path) {
        if( !path )
            throw new IllegalArgumentException("Missing 'workflow.data.store.location' configuration setting")
        if( path.startsWith(CidFileSystemProvider.SCHEME+':'))
            throw new IllegalArgumentException("Attribute 'workflow.data.store.location' cannot use 'cid:' scheme - offending value: $path")
        return FileHelper.asPath(path)
    }
}
