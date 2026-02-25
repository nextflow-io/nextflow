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

import spock.lang.Specification

/**
 * @author Edmund Miller
 */
class DatasetPathFactoryTest extends Specification {

    def 'should return null for non-dataset URIs'() {
        given:
        def factory = new DatasetPathFactory()

        expect:
        factory.parseUri('s3://bucket/key') == null
        factory.parseUri('/local/path') == null
        factory.parseUri('gs://bucket/key') == null
    }

    def 'should return null for toUriString with non-DatasetPath'() {
        given:
        def factory = new DatasetPathFactory()
        def localPath = java.nio.file.Paths.get('/tmp/test')

        expect:
        factory.toUriString(localPath) == null
    }

    def 'should return uri string for DatasetPath'() {
        given:
        def factory = new DatasetPathFactory()
        def fs = new DatasetFileSystem(new DatasetFileSystemProvider(), null)
        def path = new DatasetPath(fs, 'my-data')

        expect:
        factory.toUriString(path) == 'dataset://my-data'
    }

    def 'should return null for getBashLib'() {
        given:
        def factory = new DatasetPathFactory()
        def fs = new DatasetFileSystem(new DatasetFileSystemProvider(), null)
        def path = new DatasetPath(fs, 'my-data')

        expect:
        factory.getBashLib(path) == null
    }

    def 'should return null for getUploadCmd'() {
        given:
        def factory = new DatasetPathFactory()
        def fs = new DatasetFileSystem(new DatasetFileSystemProvider(), null)
        def path = new DatasetPath(fs, 'my-data')

        expect:
        factory.getUploadCmd('/tmp/file', path) == null
    }
}
