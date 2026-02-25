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
class DatasetPathTest extends Specification {

    DatasetFileSystem makeFs() {
        def provider = new DatasetFileSystemProvider()
        new DatasetFileSystem(provider, null)
    }

    // -- URI construction --

    def 'should parse dataset URI with host-style name'() {
        given:
        def uri = new URI('dataset://my-samplesheet')
        def path = new DatasetPath(makeFs(), uri)

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == null
        path.toString() == 'dataset://my-samplesheet'
    }

    def 'should parse dataset URI with triple-slash form'() {
        given:
        def uri = new URI('dataset:///my-samplesheet')
        def path = new DatasetPath(makeFs(), uri)

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == null
    }

    def 'should parse dataset URI with version query param'() {
        given:
        def uri = new URI('dataset://my-samplesheet?version=3')
        def path = new DatasetPath(makeFs(), uri)

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == '3'
        path.toString() == 'dataset://my-samplesheet?version=3'
    }

    // -- string construction --

    def 'should parse string path'() {
        given:
        def path = new DatasetPath(makeFs(), 'my-samplesheet')

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == null
    }

    def 'should parse string path with leading slash'() {
        given:
        def path = new DatasetPath(makeFs(), '/my-samplesheet')

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == null
    }

    def 'should parse string path with version suffix'() {
        given:
        def path = new DatasetPath(makeFs(), 'my-samplesheet@2')

        expect:
        path.datasetName == 'my-samplesheet'
        path.version == '2'
    }

    // -- Path interface --

    def 'should be absolute'() {
        expect:
        new DatasetPath(makeFs(), 'test').isAbsolute()
    }

    def 'should have name count of 1'() {
        expect:
        new DatasetPath(makeFs(), 'test').getNameCount() == 1
    }

    def 'should return self for getName(0)'() {
        given:
        def path = new DatasetPath(makeFs(), 'test')

        expect:
        path.getName(0) == path
    }

    def 'should throw for getName with invalid index'() {
        when:
        new DatasetPath(makeFs(), 'test').getName(1)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should return correct URI'() {
        given:
        def path = new DatasetPath(makeFs(), 'my-data')

        expect:
        path.toUri() == new URI('dataset://my-data')
    }

    def 'should return self for toAbsolutePath'() {
        given:
        def path = new DatasetPath(makeFs(), 'test')

        expect:
        path.toAbsolutePath().is(path)
    }

    // -- equality --

    def 'should be equal for same name and version'() {
        given:
        def a = new DatasetPath(makeFs(), 'data')
        def b = new DatasetPath(makeFs(), 'data')

        expect:
        a == b
        a.hashCode() == b.hashCode()
    }

    def 'should not be equal for different names'() {
        given:
        def a = new DatasetPath(makeFs(), 'data1')
        def b = new DatasetPath(makeFs(), 'data2')

        expect:
        a != b
    }

    def 'should not be equal for different versions'() {
        given:
        def a = new DatasetPath(makeFs(), 'data@1')
        def b = new DatasetPath(makeFs(), 'data@2')

        expect:
        a != b
    }

    // -- compareTo --

    def 'should compare by name then version'() {
        given:
        def a = new DatasetPath(makeFs(), 'alpha')
        def b = new DatasetPath(makeFs(), 'beta')
        def c = new DatasetPath(makeFs(), 'alpha@2')

        expect:
        a.compareTo(b) < 0
        b.compareTo(a) > 0
        a.compareTo(c) < 0  // null version < "2"
    }
}
