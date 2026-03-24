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
import nextflow.script.params.InParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.ValueInParam
import spock.lang.Specification

/**
 * Tests for ProcessEntryHandler parameter mapping functionality
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessEntryHandlerTest extends Specification {

    def 'should parse complex parameters with dot notation' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def result = handler.parseComplexParameters([
            'meta.id': 'SAMPLE_001',
            'meta.name': 'TestSample',
            'meta.other': 'some-value',
            'fasta': '/path/to/file.fa'
        ])

        then:
        result.meta instanceof Map
        result.meta.id == 'SAMPLE_001'
        result.meta.name == 'TestSample'
        result.meta.other == 'some-value'
        result.fasta == '/path/to/file.fa'
    }

    def 'should parse nested dot notation parameters' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def result = handler.parseComplexParameters([
            'meta.sample.id': '123',
            'meta.sample.name': 'test',
            'meta.config.quality': 'high',
            'output.dir': '/results'
        ])

        then:
        result.meta instanceof Map
        result.meta.sample instanceof Map
        result.meta.sample.id == '123'
        result.meta.sample.name == 'test'
        result.meta.config instanceof Map
        result.meta.config.quality == 'high'
        result.output instanceof Map
        result.output.dir == '/results'
    }

    def 'should handle simple parameters without dots' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def result = handler.parseComplexParameters([
            'sampleId': 'SAMPLE_001',
            'threads': '4',
            'dataFile': '/path/to/data.txt'
        ])

        then:
        result.sampleId == 'SAMPLE_001'
        result.threads == '4'
        result.dataFile == '/path/to/data.txt'
    }

    def 'should get value for val input type' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def complexParams = [
            'meta': [id: 'SAMPLE_001', name: 'TestSample'],
            'sampleId': 'SIMPLE_001'
        ]
        def valInput = Mock(ValueInParam) { getName() >> 'meta' }
        def simpleInput = Mock(ValueInParam) { getName() >> 'sampleId' }

        then:
        handler.getValueForInput(valInput, complexParams) == [id: 'SAMPLE_001', name: 'TestSample']
        handler.getValueForInput(simpleInput, complexParams) == 'SIMPLE_001'
    }

    def 'should get value for path input type' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def complexParams = [
            'fasta': '/path/to/file.fa',
            'dataFile': 'data.txt'
        ]
        def pathInput = Mock(FileInParam) { getName() >> 'fasta' }
        def fileInput = Mock(FileInParam) { getName() >> 'dataFile' }

        then:
        def fastaResult = handler.getValueForInput(pathInput, complexParams)
        def fileResult = handler.getValueForInput(fileInput, complexParams)

        // Should convert string paths to Path objects (mocked here)
        fastaResult.toString().contains('file.fa')
        fileResult.toString().contains('data.txt')
    }

    def 'should throw exception for missing required parameter' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        def complexParams = [
            'meta': [id: 'SAMPLE_001']
        ]
        def missingInput = Mock(ValueInParam) { getName() >> 'missing' }
        and:
        handler.getValueForInput(missingInput, complexParams)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should map tuple input structure correctly' () {
        given:
        def session = Mock(Session) {
            getParams() >> [
                'meta.id': 'SAMPLE_001',
                'meta.name': 'TestSample',
                'meta.other': 'some-value',
                'fasta': '/path/to/file.fa'
            ]
        }
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def processDef = Mock(ProcessDef)
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        // Mock input declaration for tuple val(meta), path(fasta)
        def tupleParam = Mock(TupleInParam) {
            getInner() >> [
                Mock(ValueInParam) { getName() >> 'meta' },
                Mock(FileInParam) { getName() >> 'fasta' }
            ]
        }

        // Test the parameter mapping logic manually
        def namedArgs = handler.parseComplexParameters(session.getParams())
        def tupleElements = []

        for( def innerParam : tupleParam.inner ) {
            def value = handler.getValueForInput(innerParam, namedArgs)
            tupleElements.add(value)
        }

        then:
        namedArgs.meta instanceof Map
        namedArgs.meta.id == 'SAMPLE_001'
        namedArgs.meta.name == 'TestSample'
        namedArgs.meta.other == 'some-value'
        namedArgs.fasta == '/path/to/file.fa'

        tupleElements.size() == 2
        tupleElements[0] instanceof Map  // meta as map
        tupleElements[0].id == 'SAMPLE_001'
        tupleElements[0].name == 'TestSample'
        tupleElements[0].other == 'some-value'
        // tupleElements[1] should be a Path object (mocked)
        tupleElements[1].toString().contains('file.fa')
    }

    def 'validateAndCastParam should pass map value through for map type' () {
        given:
        def value = [id: 'sample1', name: 'test']

        expect:
        ProcessEntryHandler.validateAndCastParam('meta', value, 'map') == value
    }

    def 'validateAndCastParam should throw for map type when value is not a map' () {
        when:
        ProcessEntryHandler.validateAndCastParam('meta', 'scalar', 'map')

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('--meta')
        e.message.contains('map')
    }

    def 'validateAndCastParam should cast string to integer' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('threads', VALUE, 'integer') == EXPECTED

        where:
        VALUE   | EXPECTED
        '4'     | 4
        '0'     | 0
        '-1'    | -1
        '99999999999' | 99999999999L
    }

    def 'validateAndCastParam should pass through existing integer value' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('threads', 4, 'integer') == 4
        ProcessEntryHandler.validateAndCastParam('threads', 4L, 'integer') == 4L
    }

    def 'validateAndCastParam should throw for integer type when value is non-numeric string' () {
        when:
        ProcessEntryHandler.validateAndCastParam('threads', 'abc', 'integer')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should throw for integer type when value is a map' () {
        when:
        ProcessEntryHandler.validateAndCastParam('threads', [a: 1], 'integer')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should cast string to float' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('ratio', '1.5', 'float') == 1.5f
        ProcessEntryHandler.validateAndCastParam('ratio', '3.14', 'float') instanceof Number
    }

    def 'validateAndCastParam should widen integer to double for float type' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('ratio', 2, 'float') == 2.0d
    }

    def 'validateAndCastParam should throw for float type when value is non-numeric string' () {
        when:
        ProcessEntryHandler.validateAndCastParam('ratio', 'abc', 'float')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should cast string to boolean' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('flag', 'true', 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', 'false', 'boolean') == Boolean.FALSE
        ProcessEntryHandler.validateAndCastParam('flag', 'TRUE', 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', 'False', 'boolean') == Boolean.FALSE
    }

    def 'validateAndCastParam should pass through existing boolean value' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('flag', Boolean.TRUE, 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', Boolean.FALSE, 'boolean') == Boolean.FALSE
    }

    def 'validateAndCastParam should throw for boolean type when value is not true or false' () {
        when:
        ProcessEntryHandler.validateAndCastParam('flag', '1', 'boolean')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should throw for boolean type when value is a map' () {
        when:
        ProcessEntryHandler.validateAndCastParam('flag', [a: 1], 'boolean')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should pass through string value for string type' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('name', 'hello', 'string') == 'hello'
    }

    def 'validateAndCastParam should throw for string type when value is a map' () {
        when:
        ProcessEntryHandler.validateAndCastParam('name', [a: 1], 'string')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should pass through path string for file and directory types' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('reads', '/path/to/file.txt', 'file') == '/path/to/file.txt'
        ProcessEntryHandler.validateAndCastParam('outdir', '/results', 'directory') == '/results'
    }

    def 'validateAndCastParam should throw for file type when value is a map' () {
        when:
        ProcessEntryHandler.validateAndCastParam('reads', [a: 1], 'file')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should pass through value unchanged for unknown type' () {
        expect:
        ProcessEntryHandler.validateAndCastParam('x', 'anything', 'unknowntype') == 'anything'
    }

    def 'should parse file input correctly' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta)
        def handler = new ProcessEntryHandler(script, session, meta)

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
}
