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

import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

/**
 * @author Edmund Miller
 */
class DatasetPathTest extends Specification {

    @Shared
    DatasetFileSystem fileSystem = new DatasetFileSystem(new DatasetFileSystemProvider(), null)

    @Unroll
    def 'should parse dataset URI #uriString'() {
        given:
        def path = new DatasetPath(fileSystem, new URI(uriString))

        expect:
        path.datasetName == expectedName
        path.version == expectedVersion
        path.toString() == expectedToString

        where:
        uriString                               | expectedName      | expectedVersion | expectedToString
        'dataset://my-samplesheet'             | 'my-samplesheet'  | null            | 'dataset://my-samplesheet'
        'dataset:///my-samplesheet'            | 'my-samplesheet'  | null            | 'dataset://my-samplesheet'
        'dataset://my-samplesheet?version=3'   | 'my-samplesheet'  | '3'             | 'dataset://my-samplesheet?version=3'
    }

    @Unroll
    def 'should parse string path #rawPath'() {
        given:
        def path = new DatasetPath(fileSystem, rawPath)

        expect:
        path.datasetName == expectedName
        path.version == expectedVersion

        where:
        rawPath             | expectedName      | expectedVersion
        'my-samplesheet'    | 'my-samplesheet'  | null
        '/my-samplesheet'   | 'my-samplesheet'  | null
        'my-samplesheet@2'  | 'my-samplesheet'  | '2'
    }

    def 'should expose basic path semantics'() {
        given:
        def path = new DatasetPath(fileSystem, 'test')

        expect:
        path.isAbsolute()
        path.getNameCount() == 1
        path.getName(0) == path
        path.toUri() == new URI('dataset://test')
        path.toAbsolutePath().is(path)
    }

    @Unroll
    def 'should throw for getName(#index)'() {
        when:
        new DatasetPath(fileSystem, 'test').getName(index)

        then:
        thrown(IllegalArgumentException)

        where:
        index << [-1, 1]
    }

    @Unroll
    def 'should compare equality for #left vs #right'() {
        given:
        def a = new DatasetPath(fileSystem, left)
        def b = new DatasetPath(fileSystem, right)

        expect:
        (a == b) == expectedEqual

        and:
        if (expectedEqual) {
            assert a.hashCode() == b.hashCode()
        }

        where:
        left      | right     | expectedEqual
        'data'    | 'data'    | true
        'data1'   | 'data2'   | false
        'data@1'  | 'data@2'  | false
    }

    @Unroll
    def 'should compare order for #left compared to #right'() {
        given:
        def a = new DatasetPath(fileSystem, left)
        def b = new DatasetPath(fileSystem, right)

        expect:
        Integer.signum(a.compareTo(b)) == expectedSign

        where:
        left      | right      | expectedSign
        'alpha'   | 'beta'     | -1
        'beta'    | 'alpha'    | 1
        'alpha'   | 'alpha@2'  | -1
    }
}
