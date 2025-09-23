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

package nextflow.script

import nextflow.Session
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
        def valInput = [type: 'val', name: 'meta']
        def simpleInput = [type: 'val', name: 'sampleId']

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
        def pathInput = [type: 'path', name: 'fasta']
        def fileInput = [type: 'file', name: 'dataFile']

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
        def missingInput = [type: 'val', name: 'missing']
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
        // Mock input structures for tuple val(meta), path(fasta)
        def inputStructures = [
            [
                type: 'tuple',
                elements: [
                    [type: 'val', name: 'meta'],
                    [type: 'path', name: 'fasta']
                ]
            ]
        ]

        // Test the parameter mapping logic manually 
        def complexParams = handler.parseComplexParameters(session.getParams())
        def tupleInput = inputStructures[0]
        def tupleElements = []
        
        for( def element : tupleInput.elements ) {
            def value = handler.getValueForInput(element, complexParams)
            tupleElements.add(value)
        }

        then:
        complexParams.meta instanceof Map
        complexParams.meta.id == 'SAMPLE_001'
        complexParams.meta.name == 'TestSample'
        complexParams.meta.other == 'some-value'
        complexParams.fasta == '/path/to/file.fa'
        
        tupleElements.size() == 2
        tupleElements[0] instanceof Map  // meta as map
        tupleElements[0].id == 'SAMPLE_001'
        tupleElements[0].name == 'TestSample'
        tupleElements[0].other == 'some-value'
        // tupleElements[1] should be a Path object (mocked)
        tupleElements[1].toString().contains('file.fa')
    }
}