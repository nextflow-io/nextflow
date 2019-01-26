/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStoragePath

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilesExTest2 extends Specification {

    @Unroll
    def 'should return uri string for #PATH' () {

        when:
        def path = PATH as Path
        then:
        path instanceof CloudStoragePath
        FilesEx.toUriString(path) == PATH

        where:
        PATH                    | _
        'gs://foo/bar'          | _
        'gs://foo'              | _
        'gs://foo/'             | _
        'gs://foo/bar/baz'      | _
        'gs://foo/bar/baz/'     | _
    }
}
