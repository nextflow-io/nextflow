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

package nextflow.script.dsl

import spock.lang.Specification

import nextflow.exception.IllegalDirectiveException
import nextflow.script.ProcessConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessBuilderTest extends Specification {

    def createBuilder() {
        new ProcessBuilder(new ProcessConfig([:]))
    }

    def 'should set directives' () {

        setup:
        def builder = createBuilder()
        def config = builder.getConfig()

        // setting list values
        when:
        builder.tag 1,2,3
        then:
        config.tag == [1,2,3]

        // setting named parameters attribute
        when:
        builder.tag field1:'val1', field2: 'val2'
        then:
        config.tag == [field1:'val1', field2: 'val2']

        // maxDuration property
        when:
        builder.time '1h'
        then:
        config.time == '1h'
        config.createTaskConfig().time == new Duration('1h')

        // maxMemory property
        when:
        builder.memory '2GB'
        then:
        config.memory == '2GB'
        config.createTaskConfig().memory == new MemoryUnit('2GB')

        when:
        builder.stageInMode 'copy'
        builder.stageOutMode 'move'
        then:
        config.stageInMode == 'copy'
        config.stageOutMode == 'move'

    }

    def 'should create PublishDir object' () {

        setup:
        def builder = createBuilder()
        def config = builder.getConfig()

        when:
        builder.publishDir '/data'
        then:
        config.get('publishDir').last() == [path:'/data']

        when:
        builder.publishDir '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir').last() == [path: '/data', mode: 'link', pattern: '*.bam']

        when:
        builder.publishDir path: '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir').last() == [path: '/data', mode: 'link', pattern: '*.bam']
    }

    def 'should throw IllegalDirectiveException'() {

        given:
        def builder = createBuilder()

        when:
        builder.hello 'world'

        then:
        def e = thrown(IllegalDirectiveException)
        e.message ==
                '''
                Unknown process directive: `hello`

                Did you mean one of these?
                        shell
                '''
                .stripIndent().trim()
    }

    def 'should set process secret'() {
        when:
        def builder = createBuilder()
        def config = builder.getConfig()
        then:
        config.getSecret() == []

        when:
        builder.secret 'foo'
        then:
        config.getSecret() == ['foo']

        when:
        builder.secret 'bar'
        then:
        config.secret == ['foo', 'bar']
        config.getSecret() == ['foo', 'bar']
    }

    def 'should set process labels'() {
        when:
        def builder = createBuilder()
        def config = builder.getConfig()
        then:
        config.getLabels() == []

        when:
        builder.label 'foo'
        then:
        config.getLabels() == ['foo']

        when:
        builder.label 'bar'
        then:
        config.getLabels() == ['foo','bar']
    }

    def 'should apply resource labels config' () {
        given:
        def builder = createBuilder()
        def config = builder.getConfig()
        expect:
        config.getResourceLabels() == [:]

        when:
        builder.resourceLabels foo: 'one', bar: 'two'
        then:
        config.getResourceLabels() == [foo: 'one', bar: 'two']

        when:
        builder.resourceLabels foo: 'new one', baz: 'three'
        then:
        config.getResourceLabels() == [foo: 'new one', bar: 'two', baz: 'three']

    }

    def 'should check a valid label' () {

        expect:
        ProcessBuilder.isValidLabel(lbl) == result

        where:
        lbl         | result
        'foo'       | true
        'foo1'      | true
        '1foo'      | false
        '_foo'      | false
        'foo1_'     | false
        'foo_1'     | true
        'foo-1'     | false
        'foo.1'     | false
        'a'         | true
        'A'         | true
        '1'         | false
        '_'         | false
        'a=b'       | true
        'a=foo'     | true
        'a=foo_1'   | true
        'a=foo_'    | false
        '_=foo'     | false
        '=a'        | false
        'a='        | false
        'a=1'       | false

    }

    def 'should apply accelerator config' () {

        given:
        def builder = createBuilder()
        def config = builder.getConfig()

        when:
        builder.accelerator 5
        then:
        config.accelerator == [limit: 5]

        when:
        builder.accelerator request: 1, limit: 5, type: 'nvida'
        then:
        config.accelerator == [request: 1, limit: 5, type: 'nvida']

        when:
        builder.accelerator 5, type: 'nvida'
        then:
        config.accelerator == [limit: 5, type: 'nvida']

        when:
        builder.accelerator 1, limit: 5
        then:
        config.accelerator == [request: 1, limit:5]

        when:
        builder.accelerator 5, request: 1
        then:
        config.accelerator == [request: 1, limit:5]
    }

    def 'should apply disk config' () {

        given:
        def builder = createBuilder()
        def config = builder.getConfig()

        when:
        builder.disk '100 GB'
        then:
        config.disk == [request: '100 GB']

        when:
        builder.disk '375 GB', type: 'local-ssd'
        then:
        config.disk == [request: '375 GB', type: 'local-ssd']

        when:
        builder.disk request: '375 GB', type: 'local-ssd'
        then:
        config.disk == [request: '375 GB', type: 'local-ssd']
    }

    def 'should apply architecture config' () {

        given:
        def builder = createBuilder()
        def config = builder.getConfig()

        when:
        builder.arch 'linux/x86_64'
        then:
        config.arch == [name: 'linux/x86_64']

        when:
        builder.arch 'linux/x86_64', target: 'zen3'
        then:
        config.arch == [name: 'linux/x86_64', target: 'zen3']

        when:
        builder.arch name: 'linux/x86_64', target: 'zen3'
        then:
        config.arch == [name: 'linux/x86_64', target: 'zen3']
    }

    def 'should throw exception for invalid error strategy' () {
        when:
        def builder = createBuilder()
        builder.errorStrategy 'abort'

        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Unknown error strategy 'abort' ― Available strategies are: terminate,finish,ignore,retry"

    }

    def 'should not throw exception for valid error strategy or closure' () {
        when:
        def builder = createBuilder()
        builder.errorStrategy 'retry'

        then:
        noExceptionThrown()

        when:
        builder = createBuilder()
        builder.errorStrategy 'terminate'

        then:
        noExceptionThrown()

        when:
        builder = createBuilder()
        builder.errorStrategy { task.exitStatus==14 ? 'retry' : 'terminate' }

        then:
        noExceptionThrown()
    }
}
