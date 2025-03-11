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

import nextflow.file.FileHelper

import java.nio.file.Path

import groovy.transform.CompileStatic
/**
 * Model data store options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataStoreOpts {

    final Path location
    final Path logLocation

    DataStoreOpts(Map opts) {
        this.location = opts.location
            ? FileHelper.toCanonicalPath(opts.location as String)
            : Path.of('.').toAbsolutePath().normalize().resolve('data')
        this.logLocation = opts.logLocation ? FileHelper.toCanonicalPath(opts.logLocation as String) : null
    }

}
