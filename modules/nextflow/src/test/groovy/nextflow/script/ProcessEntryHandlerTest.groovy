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

package nextflow.script

import java.nio.file.Path

import nextflow.Session
import nextflow.script.params.FileInParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessInputsDef
import nextflow.script.params.v2.ProcessTupleInput
import nextflow.script.types.Record
import nextflow.util.Duration
import nextflow.util.RecordMap
import spock.lang.Specification

/**
 * Tests for ProcessEntryHandler parameter mapping functionality.
 *
 * The per-type mapping of a param value to a declared input type is covered
 * by {@link ParamsHelperTest} -- the v2 tests here only assert that a typed
 * process input routes through it. The v1 tests cover the separate,
 * meta.yml-driven resolution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessEntryHandlerTest extends Specification {

    ProcessEntryHandler handler

    def setup() {
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
        handler = new ProcessEntryHandler(script, session, meta)
    }

    def 'should get value for val input type' () {
        when:
        def complexParams = [
            'meta': [id: 'SAMPLE_001', name: 'TestSample'],
            'sampleId': 'SIMPLE_001'
        ]
        def valInput = Mock(ValueInParam) { getName() >> 'meta' }
        def simpleInput = Mock(ValueInParam) { getName() >> 'sampleId' }

        then:
        handler.getValueForInputV1(valInput, complexParams, [:]) == [id: 'SAMPLE_001', name: 'TestSample']
        handler.getValueForInputV1(simpleInput, complexParams, [:]) == 'SIMPLE_001'
    }

    def 'should get value for path input type' () {
        when:
        def complexParams = [
            'fasta': '/path/to/file.fa',
            'dataFile': 'data.txt'
        ]
        def pathInput = Mock(FileInParam) { getName() >> 'fasta' }
        def fileInput = Mock(FileInParam) { getName() >> 'dataFile' }

        then:
        def fastaResult = handler.getValueForInputV1(pathInput, complexParams, [:])
        def fileResult = handler.getValueForInputV1(fileInput, complexParams, [:])

        // Should convert string paths to Path objects (mocked here)
        fastaResult.toString().contains('file.fa')
        fileResult.toString().contains('data.txt')
    }

    def 'should support destructured record input in typed process' () {
        given: 'a ProcessTupleInput representing a record(id: String, greeting: String)'
        def idInput = new ProcessInput('id', String, false)
        def greetingInput = new ProcessInput('greeting', String, false)
        def recordParam = new ProcessTupleInput([idInput, greetingInput], Record)

        and: 'a V2 config declaring that one record input'
        def inputsDef = new ProcessInputsDef()
        inputsDef.addTupleParam([idInput, greetingInput], Record)
        def processConfig = Mock(ProcessConfigV2) {
            getInputs() >> inputsDef
        }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> processConfig
        }

        when:
        def args = ProcessEntryHandler.getProcessArguments(processDef, ['id': 'abc', 'greeting': 'hello'], (Path)null)

        then: 'a single RecordMap is returned'
        args.size() == 1
        args[0] instanceof RecordMap
        args[0].id == 'abc'
        args[0].greeting == 'hello'
    }

    static class Sample implements Record {
        String id
        String text
    }

    def 'should support record type input in typed process' () {
        given: 'a process with a record type input'
        def inputsDef = new ProcessInputsDef()
        inputsDef.addParam('sample', Sample, false)
        def processConfig = Mock(ProcessConfigV2) {
            getInputs() >> inputsDef
        }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> processConfig
        }

        when:
        def args = ProcessEntryHandler.getProcessArguments(processDef, ['sample': ['id': 'a', 'text': 'hello']], (Path)null)

        then: 'the nested map is converted into a record of the declared type'
        args.size() == 1
        args[0] instanceof Record
        args[0].id == 'a'
        args[0].text == 'hello'
    }

    def 'should convert string parameter to type declared in module spec' () {
        given:
        def param = Mock(ValueInParam) { getName() >> NAME }

        expect:
        handler.getValueForInputV1(param, [(NAME): VALUE], [(NAME): TYPE]) == EXPECTED

        where:
        NAME    | TYPE    | VALUE   | EXPECTED
        'flag'  | Boolean | 'true'  | Boolean.TRUE
        'flag'  | Boolean | 'false' | Boolean.FALSE
        'flag'  | Boolean | 'TRUE'  | Boolean.TRUE
        'count' | Integer | '42'    | 42
        'count' | Integer | '0'     | 0
        'ratio' | Number  | '3.14'  | 3.14f
        // an undeclared type leaves the value unchanged
        'name'  | null    | 'hello' | 'hello'
    }

    def 'should throw error when parameter value cannot be converted to declared type' () {
        given:
        def param = Mock(ValueInParam) { getName() >> NAME }

        when:
        handler.getValueForInputV1(param, [(NAME): VALUE], [(NAME): TYPE])

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains(MESSAGE)

        where:
        NAME    | TYPE    | VALUE           | MESSAGE
        'flag'  | Boolean | 'yes'           | 'boolean (true/false)'
        'flag'  | Boolean | '1'             | 'boolean (true/false)'
        'count' | Integer | 'not-a-number'  | 'expects an integer'
        'count' | Integer | '3.14'          | 'expects an integer'
        'ratio' | Number  | 'not-a-number'  | 'floating-point'
        'meta'  | Map     | 'foo=bar'       | 'dot notation'
    }

    def 'should throw error when non-string parameter value has incompatible type' () {
        given:
        def param = Mock(ValueInParam) { getName() >> 'count' }

        when:
        handler.getValueForInputV1(param, [count: 3.14d], [count: Integer])

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('expects a Integer')
    }

    def 'should parse file input correctly' () {
        expect:
        handler.parseFileInput(INPUT) == EXPECTED_FILES

        where:
        INPUT                                                             | EXPECTED_FILES
        '/path/to/file.txt'                                               | Path.of('/path/to/file.txt')
        and:
        '/path/to/file1.txt,/path/to/file2.txt,/path/to/file3.txt'        | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        ' /path/to/file1.txt , /path/to/file2.txt , /path/to/file3.txt '  | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        '/path/to/file1.txt,,/path/to/file2.txt, ,/path/to/file3.txt'     | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        'file1.txt,file2.txt'                                             | [Path.of('file1.txt').toAbsolutePath(), Path.of('file2.txt').toAbsolutePath()]
    }

    def 'should return empty list for missing path input (v1)' () {
        given:
        def pathParam = Mock(FileInParam) { getName() >> 'intervals' }

        when: 'the input is not present in params for a path input'
        def result = handler.getValueForInputV1(pathParam, [:], [:])

        then: 'returns an empty list instead of throwing'
        result == []
    }

    def 'should reject an empty-string path input instead of treating it as not provided (v1)' () {
        given:
        def pathParam = Mock(FileInParam) { getName() >> 'proteins' }

        when: 'an empty value is supplied for a path input'
        handler.getValueForInputV1(pathParam, [proteins: ''], [:])

        then: 'it is NOT a stand-in for "not provided" -- an optional path is skipped by omitting the arg'
        def e = thrown(IllegalArgumentException)
        e.message.contains('cannot be')
    }

    def 'should throw error for missing val input (v1)' () {
        given:
        def valParam = Mock(ValueInParam) { getName() >> 'id' }

        when: 'a val input is absent'
        handler.getValueForInputV1(valParam, [:], [:])

        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Missing required parameter: --id'
    }

    def 'should resolve a typed input like a pipeline parameter (v2)' () {
        expect:
        // the full type matrix is covered by ParamsHelperTest
        handler.getValueForInputV2(new ProcessInput('x', Integer, false), [x: '42'], [x: '42']) == 42
        handler.getValueForInputV2(new ProcessInput('x', Duration, false), [x: '2h'], [x: '2h']) == Duration.of('2h')
    }

    def 'should reject a param value that cannot be converted to the declared input type (v2)' () {
        when:
        handler.getValueForInputV2(new ProcessInput('n50', type, false), [n50: value], [n50: value])

        then: 'the mismatch is reported up front instead of reaching the task'
        def e = thrown(IllegalArgumentException)
        e.message == "Parameter `--n50` with type ${type.simpleName} cannot be assigned to ${value} [String]"

        where:
        type    | value
        Integer | '3.7'      // a fractional value is not silently truncated
        Float   | 'abc'
    }

    def 'should resolve a config value as structured, not as command-line text (v2)' () {
        given:
        def decl = new ProcessInput('flag', Boolean, false)

        expect: 'a command-line value is parsed'
        handler.getValueForInputV2(decl, [flag: 'false'], [flag: 'false']) == false

        when: 'the same value comes from the config, where a Boolean is expected to be a Boolean'
        handler.getValueForInputV2(decl, [flag: 'false'], [:])

        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Parameter `--flag` with type Boolean cannot be assigned to false [String]'
    }

    def 'should return null for missing optional input (v2)' () {
        // Regression test for https://github.com/nextflow-io/nextflow/issues/7161
        given:
        def param = new ProcessInput('gzi', Path, true)

        when: 'a v2 optional path input is not provided'
        def result = handler.getValueForInputV2(param, [:], [:])

        then: 'returns null instead of throwing'
        result == null
    }

    def 'should throw error for missing required input (v2)' () {
        given:
        def param = new ProcessInput('reads', Path, false)

        when:
        handler.getValueForInputV2(param, [:], [:])

        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Missing required parameter: --reads'
    }
}
