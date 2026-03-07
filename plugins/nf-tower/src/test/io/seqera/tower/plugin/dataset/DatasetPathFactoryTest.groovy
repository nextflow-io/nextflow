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

package io.seqera.tower.plugin.dataset

import java.nio.file.Paths

import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

/**
 * @author Edmund Miller
 */
class DatasetPathFactoryTest extends Specification {

    @Shared
    DatasetPathFactory factory = new DatasetPathFactory()

    @Shared
    DatasetFileSystem datasetFs = new DatasetFileSystem(new DatasetFileSystemProvider(), null)

    @Unroll
    def 'parseUri should return null for non-dataset URI #value'() {
        expect:
        factory.parseUri(value) == null

        where:
        value << ['s3://bucket/key', '/local/path', 'gs://bucket/key']
    }

    @Unroll
    def 'toUriString should return #expected for #scenario'() {
        given:
        def path = pathFactory.call()

        expect:
        factory.toUriString(path) == expected

        where:
        scenario              | pathFactory                                         | expected
        'DatasetPath input'   | { new DatasetPath(datasetFs, 'my-data') }          | 'dataset://my-data'
        'non-DatasetPath'     | { Paths.get('/tmp/test') }                         | null
    }

    @Unroll
    def '#methodName should return null for dataset path'() {
        given:
        def path = new DatasetPath(datasetFs, 'my-data')

        expect:
        invoke.call(factory, path) == null

        where:
        methodName     | invoke
        'getBashLib'   | { DatasetPathFactory f, target -> f.getBashLib(target) }
        'getUploadCmd' | { DatasetPathFactory f, target -> f.getUploadCmd('/tmp/file', target) }
    }
}
