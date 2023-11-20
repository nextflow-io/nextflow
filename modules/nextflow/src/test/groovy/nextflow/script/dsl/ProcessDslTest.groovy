/*
 * Copyright 2013-2023, Seqera Labs
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
import spock.lang.Unroll

import nextflow.exception.IllegalDirectiveException
import nextflow.script.params.FileInParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.script.TokenVar
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDslTest extends Specification {

    def 'should set directives' () {

        setup:
        def builder = new ProcessDsl(Mock(BaseScript), null)
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


    def 'should create input directives' () {

        setup:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        def config = builder.getConfig()

        when:
        builder._in_file([infile:'filename.fa'])
        builder._in_val('x').setFrom(1)
        builder._in_stdin()

        then:
        config.getInputs().size() == 3

        config.inputs.get(0) instanceof FileInParam
        config.inputs.get(0).name == 'infile'
        (config.inputs.get(0) as FileInParam).filePattern == 'filename.fa'

        config.inputs.get(1) instanceof ValueInParam
        config.inputs.get(1).name == 'x'

        config.inputs.get(2).name == '-'
        config.inputs.get(2) instanceof StdInParam

        config.inputs.names == [ 'infile', 'x', '-' ]
        config.inputs.ofType( FileInParam ) == [ config.getInputs().get(0) ]

    }

    def 'should create output directives' () {

        setup:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        def config = builder.getConfig()

        when:
        builder._out_stdout()
        builder._out_file(new TokenVar('file1')).setInto('ch1')
        builder._out_file(new TokenVar('file2')).setInto('ch2')
        builder._out_file(new TokenVar('file3')).setInto('ch3')

        then:
        config.outputs.size() == 4
        config.outputs.names == ['-', 'file1', 'file2', 'file3']
        config.outputs.ofType(StdOutParam).size() == 1

        config.outputs[0] instanceof StdOutParam
        config.outputs[1].name == 'file1'
        config.outputs[2].name == 'file2'
        config.outputs[3].name == 'file3'

    }

    def 'should create PublishDir object' () {

        setup:
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        def builder = new ProcessDsl(Mock(BaseScript), null)

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
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        ProcessDsl.isValidLabel(lbl) == result

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

    @Unroll
    def 'should match selector: #SELECTOR with #TARGET' () {
        expect:
        ProcessDsl.matchesSelector(TARGET, SELECTOR) == EXPECTED

        where:
        SELECTOR        | TARGET    | EXPECTED
        'foo'           | 'foo'     | true
        'foo'           | 'bar'     | false
        '!foo'          | 'bar'     | true
        'a|b'           | 'a'       | true
        'a|b'           | 'b'       | true
        'a|b'           | 'z'       | false
        'a*'            | 'a'       | true
        'a*'            | 'aaaa'    | true
        'a*'            | 'bbbb'    | false
    }

    def 'should apply config setting for a process label' () {
        given:
        def settings = [
                'withLabel:short'  : [ cpus: 1, time: '1h'],
                'withLabel:!short' : [ cpus: 32, queue: 'cn-long'],
                'withLabel:foo'    : [ cpus: 2 ],
                'withLabel:foo|bar': [ disk: '100GB' ],
                'withLabel:gpu.+'  : [ cpus: 4 ],
        ]

        when:
        def config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithLabels(settings, ['short'])
        then:
        config.cpus == 1
        config.time == '1h'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithLabels(settings, ['long'])
        then:
        config.cpus == 32
        config.queue == 'cn-long'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithLabels(settings, ['foo'])
        then:
        config.cpus == 2
        config.disk == '100GB'
        config.queue == 'cn-long'
        config.size() == 3

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithLabels(settings, ['bar'])
        then:
        config.cpus == 32
        config.disk == '100GB'
        config.queue == 'cn-long'
        config.size() == 3

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithLabels(settings, ['gpu-1'])
        then:
        config.cpus == 4
        config.queue == 'cn-long'
        config.size() == 2

    }


    def 'should apply config setting for a process name' () {
        given:
        def settings = [
                'withName:alpha'        : [ cpus: 1, time: '1h'],
                'withName:delta'        : [ cpus: 2 ],
                'withName:delta|gamma'  : [ disk: '100GB' ],
                'withName:omega.+'      : [ cpus: 4 ],
        ]

        when:
        def config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithName(settings, 'xx')
        then:
        config.size() == 0

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithName(settings, 'alpha')
        then:
        config.cpus == 1
        config.time == '1h'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithName(settings, 'delta')
        then:
        config.cpus == 2
        config.disk == '100GB'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithName(settings, 'gamma')
        then:
        config.disk == '100GB'
        config.size() == 1

        when:
        config = new ProcessConfig([:])
        new ProcessDsl(config).applyConfigSelectorWithName(settings, 'omega_x')
        then:
        config.cpus == 4
        config.size() == 1
    }


    def 'should apply config process defaults' () {

        when:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        builder.queue 'cn-el6'
        builder.memory '10 GB'
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest'
        )
        def config = builder.getConfig()

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.cacheable == true



        when:
        builder = new ProcessDsl(Mock(BaseScript), null)
        builder.container null
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        config = builder.getConfig()
        then:
        config.queue == 'def-queue'
        config.container == null
        config.maxRetries == 5



        when:
        builder = new ProcessDsl(Mock(BaseScript), null)
        builder.maxRetries 10
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        config = builder.getConfig()
        then:
        config.queue == 'def-queue'
        config.container == 'ubuntu:latest'
        config.maxRetries == 10
    }


    def 'should apply pod configs' () {

        when:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        builder.applyConfigDefaults( pod: [secret: 'foo', mountPath: '/there'] )
        then:
        builder.getConfig().pod == [
                [secret: 'foo', mountPath: '/there']
        ]

        when:
        builder = new ProcessDsl(Mock(BaseScript), null)
        builder.applyConfigDefaults( pod: [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ] )
        then:
        builder.getConfig().pod == [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ]

    }

    def 'should clone config object' () {

        when:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        def config = builder.getConfig()
        builder.queue 'cn-el6'
        builder.container 'ubuntu:latest'
        builder.memory '10 GB'
        builder._in_val('foo')
        builder._in_file('sample.txt')
        builder._out_file('result.txt')

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().size() == 2
        config.getOutputs().size() == 1

        when:
        def copy = config.clone()
        builder = new ProcessDsl(copy)
        builder.queue 'long'
        builder.container 'debian:wheezy'
        builder.memory '5 GB'
        builder._in_val('bar')
        builder._out_file('sample.bam')

        then:
        copy.queue == 'long'
        copy.container == 'debian:wheezy'
        copy.memory == '5 GB'
        copy.getInputs().size() == 3
        copy.getOutputs().size() == 2

        // original config is not affected
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().size() == 2
        config.getOutputs().size() == 1
    }

    def 'should apply accelerator config' () {

        given:
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        def builder = new ProcessDsl(Mock(BaseScript), null)
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
        def builder = new ProcessDsl(Mock(BaseScript), null)
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

    def 'should not apply config on negative label' () {
        given:
        def settings = [
                'withLabel:foo': [ cpus: 2 ],
                'withLabel:!foo': [ cpus: 4 ],
                'withLabel:!nodisk_.*': [ disk: '100.GB']
        ]

        when:
        def config = new ProcessConfig(label: ['foo', 'other'])
        new ProcessDsl(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 2
        config.disk == '100.GB'

        when:
        config = new ProcessConfig(label: ['foo', 'other', 'nodisk_label'])
        new ProcessDsl(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 2
        !config.disk

        when:
        config = new ProcessConfig(label: ['other', 'nodisk_label'])
        new ProcessDsl(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 4
        !config.disk

    }

    def 'should throw exception for invalid error strategy' () {
        when:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        builder.errorStrategy 'abort'

        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Unknown error strategy 'abort' ― Available strategies are: terminate,finish,ignore,retry"

    }

    def 'should not throw exception for valid error strategy or closure' () {
        when:
        def builder = new ProcessDsl(Mock(BaseScript), null)
        builder.errorStrategy 'retry'

        then:
        noExceptionThrown()

        when:
        builder = new ProcessDsl(Mock(BaseScript), null)
        builder.errorStrategy 'terminate'

        then:
        noExceptionThrown()

        when:
        builder = new ProcessDsl(Mock(BaseScript), null)
        builder.errorStrategy { task.exitStatus==14 ? 'retry' : 'terminate' }

        then:
        noExceptionThrown()
    }
}
