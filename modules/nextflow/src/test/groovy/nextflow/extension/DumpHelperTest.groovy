/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.yaml.YamlSlurper
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class DumpHelperTest extends Specification {

    def 'should normalize a value'() {
        given:
        def path = Path.of('/some/file.txt')

        expect: 'normalize null values to the given sentinel'
        DumpHelper.normalize(null) == null
        DumpHelper.normalize(null, 'N/A') == 'N/A'
        DumpHelper.normalize(null, 0) == 0

        and:
        DumpHelper.normalize(42) == 42
        DumpHelper.normalize(3.14) == 3.14
        DumpHelper.normalize(true) == true
        DumpHelper.normalize('hello') == 'hello'

        and: 'normalize a path to its URI string'
        DumpHelper.normalize(path) == path.toString()

        and: 'normalize a collection recursively'
        DumpHelper.normalize([1, 'two', null, path]) == [1, 'two', null, path.toString()]
        DumpHelper.normalize([[1, 2], [3, 4]]) == [[1, 2], [3, 4]]

        and: 'normalize a map recursively'
        DumpHelper.normalize([a: 1, b: path, c: null]) == [a: 1, b: path.toString(), c: null]
        DumpHelper.normalize([outer: [inner: 42]]) == [outer: [inner: 42]]
        DumpHelper.normalize([files: [path]]) == [files: [path.toString()]]
    }

    def 'should pretty-print a value'() {
        when:
        def value = ['a', 'b', 'c']
        def result = DumpHelper.prettyPrint(value)
        then:
        result.contains('"a"')
        result.contains('"b"')
        result.contains('"c"')

        when:
        value = [key: 'value']
        result = DumpHelper.prettyPrint(value)
        then:
        result.contains('"key"')
        result.contains('"value"')

        when:
        value = Path.of('/some/path')
        result = DumpHelper.prettyPrint(value)
        then:
        result.contains(value.toString())

        expect:
        DumpHelper.prettyPrint('hello') == "'hello'"
        DumpHelper.prettyPrint(42) == '42'
        DumpHelper.prettyPrint(true) == 'true'
    }

    def 'should produce pretty-printed JSON'() {

        when: 'for a collection'
        def value = [1, 2, 3]
        def result = DumpHelper.prettyPrintJson(value)
        then:
        fromJson(result) == value

        when: 'for a map'
        value = [name: 'Alice', age: 30]
        result = DumpHelper.prettyPrintJson(value)
        then:
        fromJson(result) == value

        when: 'for a path'
        value = Path.of('/data/sample.bam')
        result = DumpHelper.prettyPrintJson([key: value])
        then:
        fromJson(result) == [key: value.toString()]

        when: 'for null values'
        value = [key: null]
        result = DumpHelper.prettyPrintJson(value)
        then:
        fromJson(result) == value
    }

    def 'should produce pretty-printed YAML'() {
        when: 'for a collection'
        def value = ['a', 'b', 'c']
        def result = DumpHelper.prettyPrintYaml(value)
        then:
        fromYaml(result) == value

        when: 'for a map'
        value = [name: 'Alice', age: 30]
        result = DumpHelper.prettyPrintYaml(value)
        then:
        fromYaml(result) == value

        when: 'for a path'
        value = Path.of('/data/sample.bam')
        result = DumpHelper.prettyPrintYaml([file: value])
        then:
        fromYaml(result) == [file: value.toString()]
    }

    private Object fromJson(String text) {
        return new JsonSlurper().parseText(text)
    }

    private Object fromYaml(String text) {
        return new YamlSlurper().parseText(text)
    }
}
