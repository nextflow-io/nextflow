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

import nextflow.script.params.InputsList

import java.nio.file.Files
import java.nio.file.Path
import nextflow.Session
import nextflow.script.params.FileInParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessInputsDef
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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

    def 'should handle dotted key overwriting plain value' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
        def handler = new ProcessEntryHandler(script, session, meta)

        when:
        'both --foo bar and --foo.id 1 are passed, the dotted key wins'
        def result = handler.parseComplexParameters([
            'foo': 'bar',
            'foo.id': '1'
        ])

        then:
        result.foo instanceof Map
        result.foo.id == '1'
    }

    def 'should get value for val input type' () {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> [ 'hello' ]
        }
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

    def 'validateAndCastParam should pass map value through for map type'() {
        given:
        def value = [id: 'sample1', name: 'test']

        expect:
        ProcessEntryHandler.validateAndCastParam('meta', value, 'map') == value
    }

    def 'validateAndCastParam should throw for map type when value is not a map'() {
        when:
        ProcessEntryHandler.validateAndCastParam('meta', 'scalar', 'map')

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('--meta')
        e.message.contains('map')
    }

    def 'validateAndCastParam should cast string to integer'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('threads', VALUE, 'integer') == EXPECTED

        where:
        VALUE         | EXPECTED
        '4'           | 4
        '0'           | 0
        '-1'          | -1
        '99999999999' | 99999999999L
    }

    def 'validateAndCastParam should pass through existing integer value'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('threads', 4, 'integer') == 4
        ProcessEntryHandler.validateAndCastParam('threads', 4L, 'integer') == 4L
    }

    def 'validateAndCastParam should throw for integer type when value is non-numeric string'() {
        when:
        ProcessEntryHandler.validateAndCastParam('threads', 'abc', 'integer')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should throw for integer type when value is a map'() {
        when:
        ProcessEntryHandler.validateAndCastParam('threads', [a: 1], 'integer')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should cast string to float'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('ratio', '1.5', 'float') == 1.5f
        ProcessEntryHandler.validateAndCastParam('ratio', '3.14', 'float') instanceof Number
    }

    def 'validateAndCastParam should widen integer to double for float type'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('ratio', 2, 'float') == 2.0d
    }

    def 'validateAndCastParam should throw for float type when value is non-numeric string'() {
        when:
        ProcessEntryHandler.validateAndCastParam('ratio', 'abc', 'float')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should cast string to boolean'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('flag', 'true', 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', 'false', 'boolean') == Boolean.FALSE
        ProcessEntryHandler.validateAndCastParam('flag', 'TRUE', 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', 'False', 'boolean') == Boolean.FALSE
    }

    def 'validateAndCastParam should pass through existing boolean value'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('flag', Boolean.TRUE, 'boolean') == Boolean.TRUE
        ProcessEntryHandler.validateAndCastParam('flag', Boolean.FALSE, 'boolean') == Boolean.FALSE
    }

    def 'validateAndCastParam should throw for boolean type when value is not true or false'() {
        when:
        ProcessEntryHandler.validateAndCastParam('flag', '1', 'boolean')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should throw for boolean type when value is a map'() {
        when:
        ProcessEntryHandler.validateAndCastParam('flag', [a: 1], 'boolean')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should pass through string value for string type'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('name', 'hello', 'string') == 'hello'
    }

    def 'validateAndCastParam should pass through path string for file and directory types'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('reads', '/path/to/file.txt', 'file') == Path.of('/path/to/file.txt')
        ProcessEntryHandler.validateAndCastParam('outdir', '/results', 'directory') == Path.of('/results')
    }

    def 'validateAndCastParam should throw for file type when value is a map'() {
        when:
        ProcessEntryHandler.validateAndCastParam('reads', [a: 1], 'file')

        then:
        thrown(IllegalArgumentException)
    }

    def 'validateAndCastParam should pass through value unchanged for unknown type'() {
        expect:
        ProcessEntryHandler.validateAndCastParam('x', 'anything', 'unknowntype') == 'anything'
    }

    def 'should parse file input correctly'() {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def meta = Mock(ScriptMeta) {
            getLocalProcessNames() >> ['hello']
        }
        def handler = new ProcessEntryHandler(script, session, meta)

        expect:
        handler.parseFileInput(INPUT) == EXPECTED_FILES

        where:
        INPUT                                                            | EXPECTED_FILES
        '/path/to/file.txt'                                              | Path.of('/path/to/file.txt')
        and:
        '/path/to/file1.txt,/path/to/file2.txt,/path/to/file3.txt'       | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        ' /path/to/file1.txt , /path/to/file2.txt , /path/to/file3.txt ' | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        '/path/to/file1.txt,,/path/to/file2.txt, ,/path/to/file3.txt'    | [Path.of('/path/to/file1.txt'), Path.of('/path/to/file2.txt'), Path.of('/path/to/file3.txt')]
        'file1.txt,file2.txt'                                            | [Path.of('file1.txt').toAbsolutePath(), Path.of('file2.txt').toAbsolutePath()]
    }

    // --- getProcessArguments: V1 config ---

    def 'should get process arguments for V1 config with simple val inputs'() {
        given:
        def session = Mock(Session) {
            getParams() >> [sampleId: 'S001', threads: '4']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.addAll([
            Mock(ValueInParam) { getName() >> 'sampleId' },
            Mock(ValueInParam) { getName() >> 'threads' }
        ])
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result == ['S001', '4']
    }

    def 'should get process arguments for V1 config with file input'() {
        given:
        def session = Mock(Session) {
            getParams() >> [reads: '/data/sample.fastq']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.add(Mock(FileInParam) { getName() >> 'reads' })
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof Path
        result[0].toString().endsWith('sample.fastq')
    }

    def 'should get process arguments for V1 config with comma-separated file list'() {
        given:
        def session = Mock(Session) {
            getParams() >> [reads: '/data/r1.fastq,/data/r2.fastq']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.add(Mock(FileInParam) { getName() >> 'reads' })
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof List
        def files = result[0] as List
        files.size() == 2
        files[0] instanceof Path
        files[0].toString().endsWith('r1.fastq')
        files[1] instanceof Path
        files[1].toString().endsWith('r2.fastq')
    }

    def 'should get process arguments for V1 config with map parameter from dot notation'() {
        given:
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', 'meta.name': 'Test']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.add(Mock(ValueInParam) { getName() >> 'meta' })

        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof Map
        (result[0] as Map).id == 'S001'
        (result[0] as Map).name == 'Test'
    }

    def 'should get process arguments for V1 config with tuple input (val + path)'() {
        given:
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', 'meta.name': 'Test', fasta: '/data/genome.fa']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def tupleParam = Mock(TupleInParam) {
            getInner() >> [
                Mock(ValueInParam) { getName() >> 'meta' },
                Mock(FileInParam) { getName() >> 'fasta' }
            ]
        }
        def inputs = new InputsList()
        inputs.add(tupleParam)
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof List
        def tuple = result[0] as List
        tuple.size() == 2
        tuple[0] instanceof Map
        (tuple[0] as Map).id == 'S001'
        (tuple[0] as Map).name == 'Test'
        tuple[1] instanceof Path
        tuple[1].toString().endsWith('genome.fa')
    }

    def 'should return empty list for V1 config with no inputs'() {
        given:
        def session = Mock(Session)
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def config = Mock(ProcessConfigV1) { getInputs() >> new InputsList() }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, [:])

        then:
        result == []
    }

    // --- getProcessArguments: V2 config ---

    private static ProcessInputsDef v2Inputs(ProcessInput... inputs) {
        def def_ = new ProcessInputsDef()
        inputs.each { def_.params.add(it) }
        return def_
    }

    def 'should get process arguments for V2 config with String input'() {
        given:
        def session = Mock(Session) {
            getParams() >> [sampleId: 'S001']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('sampleId', String, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result == ['S001']
    }

    def 'should get process arguments for V2 config with integer type casting from string'() {
        given:
        def session = Mock(Session) {
            getParams() >> [threads: '8']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('threads', Integer, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result == [8]
        result[0] instanceof Integer
    }

    def 'should get process arguments for V2 config with float type casting from string'() {
        given:
        def session = Mock(Session) {
            getParams() >> [ratio: '0.75']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('ratio', Float, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof Number
        (result[0] as Number).floatValue() == 0.75f
    }

    def 'should get process arguments for V2 config with boolean type casting from string'() {
        given:
        def session = Mock(Session) {
            getParams() >> [flag: 'true']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('flag', Boolean, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result == [Boolean.TRUE]
    }

    def 'should get process arguments for V2 config with Path input'() {
        given:
        def session = Mock(Session) {
            getParams() >> [reads: '/data/sample.fastq']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('reads', Path, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof Path
        result[0].toString().endsWith('sample.fastq')
    }

    def 'should get process arguments for V2 config with comma-separated file list'() {
        given:
        def session = Mock(Session) {
            getParams() >> [reads: '/data/r1.fastq,/data/r2.fastq,/data/r3.fastq']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('reads', List<Path>, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof List
        def files = result[0] as List
        files.size() == 3
        files[0].toString().endsWith('r1.fastq')
        files[1].toString().endsWith('r2.fastq')
        files[2].toString().endsWith('r3.fastq')
    }

    def 'should get process arguments for V2 config with map from dot notation'() {
        given:
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', 'meta.name': 'Test', 'meta.paired': 'true']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(new ProcessInput('meta', Map, false))
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof Map
        (result[0] as Map).id == 'S001'
        (result[0] as Map).name == 'Test'
        (result[0] as Map).paired == 'true'
    }

    def 'should get process arguments for V2 config with tuple (val map + path)'() {
        given:
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', 'meta.name': 'Test', fasta: '/data/genome.fa']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def tupleComponents = [
            new ProcessInput('meta', Map, false),
            new ProcessInput('fasta', Path, false)
        ]
        def inputsDef = new ProcessInputsDef()
        inputsDef.addTupleParam(tupleComponents, List)
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof List
        def tuple = result[0] as List
        tuple.size() == 2
        tuple[0] instanceof Map
        (tuple[0] as Map).id == 'S001'
        (tuple[0] as Map).name == 'Test'
        tuple[1] instanceof Path
        tuple[1].toString().endsWith('genome.fa')
    }

    def 'should get process arguments for V2 config with multiple mixed inputs'() {
        given:
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', reads: '/data/r1.fastq', threads: '4']
        }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(
            new ProcessInput('meta', Map, false),
            new ProcessInput('reads', Path, false),
            new ProcessInput('threads', Integer, false)
        )
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 3
        (result[0] as Map).id == 'S001'
        result[1] instanceof Path
        result[1].toString().endsWith('r1.fastq')
        result[2] == 4
    }

    def 'should return empty list for V2 config with no inputs'() {
        given:
        def session = Mock(Session) { getParams() >> [:] }
        def handler = new ProcessEntryHandler(Mock(BaseScript), session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def config = Mock(ProcessConfigV2) { getInputs() >> new ProcessInputsDef() }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result == []
    }

    // --- getProcessArguments with meta.yml type casting ---

    /**
     * Creates a mock BaseScript whose binding points to a fake main.nf in a temp dir
     * that also contains meta.yml with the given YAML content.
     */
    private BaseScript scriptWithMetaYml(String yamlContent) {
        def tmpDir = Files.createTempDirectory('nxf-test-meta-')
        tmpDir.toFile().deleteOnExit()
        Files.write(tmpDir.resolve('meta.yml'), yamlContent.bytes)
        def scriptPath = tmpDir.resolve('main.nf')
        Files.write(scriptPath, ''.bytes)
        def binding = Mock(ScriptBinding) { getScriptPath() >> scriptPath }
        return Mock(BaseScript) { getBinding() >> binding }
    }

    def 'should cast integer, float and boolean params via meta.yml in V1 config'() {
        given:
        def script = scriptWithMetaYml('''\
input:
  - name: sampleId
    type: string
  - name: threads
    type: integer
  - name: ratio
    type: float
  - name: paired
    type: boolean
''')
        def session = Mock(Session) {
            getParams() >> [sampleId: 'S001', threads: '8', ratio: '0.75', 'paired': 'false']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.addAll([
            Mock(ValueInParam) { getName() >> 'sampleId' },
            Mock(ValueInParam) { getName() >> 'threads' },
            Mock(ValueInParam) { getName() >> 'ratio' },
            Mock(ValueInParam) { getName() >> 'paired' }
        ])
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result[0] == 'S001'
        result[1] == 8
        result[1] instanceof Integer
            result[2] instanceof Number
        (result[2] as Number).floatValue() == 0.75f
        result[3] == Boolean.FALSE

    }

    def 'should accept meta map and file path via meta.yml tuple declaration in V1 config'() {
        given:
        // Tuple in meta.yml: list-of-list (paramSpec format)
        def script = scriptWithMetaYml('''\
input:
  - - name: meta
      type: map
    - name: fasta
      type: file
''')
        def session = Mock(Session) {
            getParams() >> ['meta.id': 'S001', 'meta.name': 'Test', fasta: '/data/genome.fa']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def tupleParam = Mock(TupleInParam) {
            getInner() >> [
                Mock(ValueInParam) { getName() >> 'meta' },
                Mock(FileInParam) { getName() >> 'fasta' }
            ]
        }
        def inputs = new InputsList()
        inputs.add(tupleParam)
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        def tuple = result[0] as List
        tuple[0] instanceof Map
        (tuple[0] as Map).id == 'S001'
        (tuple[0] as Map).name == 'Test'
        tuple[1].toString().endsWith('genome.fa')
    }

    def 'should handle comma-separated file list with meta.yml file type in V1 config'() {
        given:
        def script = scriptWithMetaYml('''\
input:
  - name: reads
    type: file
''')
        def session = Mock(Session) {
            getParams() >> [reads: '/data/r1.fastq,/data/r2.fastq,/data/r3.fastq']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.add(Mock(FileInParam) { getName() >> 'reads' })
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result.size() == 1
        result[0] instanceof List<Path>
        def files = result[0] as List
        files.size() == 3
        files[0].toString().endsWith('r1.fastq')
        files[1].toString().endsWith('r2.fastq')
        files[2].toString().endsWith('r3.fastq')
    }

    def 'should cast params via old nf-core meta.yml format in V1 config'() {
        given:
        def script = scriptWithMetaYml('''\
input:
  - sampleId:
      type: string
  - threads:
      type: integer
  - ratio:
      type: float
''')
        def session = Mock(Session) {
            getParams() >> [sampleId: 'S001', threads: '16', ratio: '1.5']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.addAll([
            Mock(ValueInParam) { getName() >> 'sampleId' },
            Mock(ValueInParam) { getName() >> 'threads' },
            Mock(ValueInParam) { getName() >> 'ratio' }
        ])
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result[0] == 'S001'
        result[1] == 16
        result[1] instanceof Integer
        (result[2] as Number).floatValue() == 1.5f
    }

    def 'should cast params via meta.yml in V2 config with typed ProcessInput'() {
        given:
        // meta.yml declares integer for threads; ProcessInput also says Integer — both agree
        def script = scriptWithMetaYml('''\
input:
  - name: sampleId
    type: string
  - name: threads
    type: integer
''')
        def session = Mock(Session) {
            getParams() >> [sampleId: 'S001', threads: '4']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputsDef = v2Inputs(
            new ProcessInput('sampleId', String, false),
            new ProcessInput('threads', Integer, false)
        )
        def config = Mock(ProcessConfigV2) { getInputs() >> inputsDef }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        def result = handler.getProcessArguments(processDef, session.params)

        then:
        result[0] == 'S001'
        result[1] == 4
        result[1] instanceof Integer
    }

    def 'should throw when meta.yml type conflicts with provided value in V1 config'() {
        given:
        def script = scriptWithMetaYml('''\
input:
  - name: threads
    type: integer
''')
        def session = Mock(Session) {
            getParams() >> [threads: 'not-a-number']
        }
        def handler = new ProcessEntryHandler(script, session, Mock(ScriptMeta) { getLocalProcessNames() >> ['hello'] })

        def inputs = new InputsList()
        inputs.add(Mock(ValueInParam) { getName() >> 'threads' })
        def config = Mock(ProcessConfigV1) { getInputs() >> inputs }
        def processDef = Mock(ProcessDef) {
            getProcessConfig() >> config
            getName() >> 'myProcess'
        }

        when:
        handler.getProcessArguments(processDef, session.params)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('--threads')
    }
}
