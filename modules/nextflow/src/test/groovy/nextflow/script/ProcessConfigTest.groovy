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

import java.nio.file.Files

import nextflow.scm.ProviderConfig
import spock.lang.Specification
import spock.lang.Unroll

import nextflow.exception.IllegalDirectiveException
import nextflow.processor.ErrorStrategy
import nextflow.script.params.FileInParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueInParam
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import static nextflow.util.CacheHelper.HashMode
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessConfigTest extends Specification {


    def 'should return defaults' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        expect:
        config.shell ==  ['/bin/bash','-ue']
        config.cacheable
        config.maxRetries == 1
        config.maxErrors == -1
        config.errorStrategy == ErrorStrategy.TERMINATE
    }

    def 'should set properties' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        // setting property using method without brackets
        when:
        config.tag = 'val 1'
        then:
        config.tag == 'val 1'

        // setting list values
        when:
        config.tag 1,2,3
        then:
        config.tag == [1,2,3]

        // setting named parameters attribute
        when:
        config.tag field1:'val1', field2: 'val2'
        then:
        config.tag == [field1:'val1', field2: 'val2']

        // generic value assigned like a 'plain' property
        when:
        config.tag = 99
        then:
        config.tag == 99

        // maxDuration property
        when:
        config.time '1h'
        then:
        config.time == '1h'
        config.createTaskConfig().time == new Duration('1h')

        // maxMemory property
        when:
        config.memory '2GB'
        then:
        config.memory == '2GB'
        config.createTaskConfig().memory == new MemoryUnit('2GB')

        when:
        config.stageInMode 'copy'
        config.stageOutMode 'move'
        then:
        config.stageInMode == 'copy'
        config.stageOutMode == 'move'

    }

    @Unroll
    def 'should set fair directive' () {
        given:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config.fair = CONFIG
        then:
        config.getFair() == EXPECTED

        where:
        CONFIG      | EXPECTED
        null        | false
        false       | false
        true        | true
    }

    def 'should parse properties'() {

        when:
        def config = new ProcessConfig( maxDuration:'1h' )
        then:
        config.maxDuration as Duration == Duration.of('1h')
    }


    def 'should not throw MissingPropertyException' () {

        when:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)
        def x = config.hola

        then:
        x == null
        noExceptionThrown()

    }

    def 'should throw MissingPropertyException' () {
        when:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script).throwExceptionOnMissingProperty(true)
        def x = config.hola

        then:
        thrown(MissingPropertyException)
    }


    def 'should check property existence' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        expect:
        config.containsKey('debug')
        config.containsKey('shell')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null

    }

    def 'should create input directives' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config._in_file([infile:'filename.fa'])
        config._in_val('x').setFrom(1)
        config._in_stdin()

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
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config._out_stdout()
        config._out_file(new TokenVar('file1')).setInto('ch1')
        config._out_file(new TokenVar('file2')).setInto('ch2')
        config._out_file(new TokenVar('file3')).setInto('ch3')

        then:
        config.outputs.size() == 4
        config.outputs.names == ['-', 'file1', 'file2', 'file3']
        config.outputs.ofType(StdOutParam).size() == 1

        config.outputs[0] instanceof StdOutParam
        config.outputs[1].name == 'file1'
        config.outputs[2].name == 'file2'
        config.outputs[3].name == 'file3'

    }


    def 'should set cache attribute'() {

        when:
        def config = new ProcessConfig(map)
        then:
        config.cacheable == result
        config.isCacheable() == result
        config.getHashMode() == mode

        where:
        result | mode               | map
        true   | HashMode.STANDARD  | [:]
        true   | HashMode.STANDARD  | [cache:true]
        true   | HashMode.STANDARD  | [cache:'yes']
        true   | HashMode.DEEP      | [cache:'deep']
        false  | HashMode.STANDARD  | [cache:false]
        false  | HashMode.STANDARD  | [cache:'false']
        false  | HashMode.STANDARD  | [cache:'off']
        false  | HashMode.STANDARD  | [cache:'no']

    }


    def 'should set ext property' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        // setting property using method without brackets
        when:
        config.ext.tool = 'blast'
        config.ext.modules = ['a/1', 'b/2']
        config.ext.command = { "echo $foo" }
        then:
        config.ext.tool == 'blast'
        config.ext.modules ==  ['a/1', 'b/2']

        when:
        def task = config.createTaskConfig()
        task.setContext( [foo: 'Hello'] )
        then:
        task.ext.tool == 'blast'
        task.ext.modules.join(',') == 'a/1,b/2'
        task.ext.command == 'echo Hello'

    }

    def 'should create PublishDir object' () {

        setup:
        BaseScript script = Mock(BaseScript)
        ProcessConfig config

        when:
        config = new ProcessConfig(script)
        config.publishDir '/data'
        then:
        config.get('publishDir')[0] == [path:'/data']

        when:
        config = new ProcessConfig(script)
        config.publishDir '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir')[0] == [path: '/data', mode: 'link', pattern: '*.bam']

        when:
        config = new ProcessConfig(script)
        config.publishDir path: '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir')[0] == [path: '/data', mode: 'link', pattern: '*.bam']
    }

    def 'should throw InvalidDirectiveException'() {

        given:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config.hello 'world'

        then:
        def e = thrown(IllegalDirectiveException)
        e.message ==
                '''
                Unknown process directive: `hello`

                Did you mean of these?
                        shell
                '''
                .stripIndent().trim()
    }

    def 'should set process secret'() {
        when:
        def config = new ProcessConfig([:])
        then:
        config.getSecret() == []

        when:
        config.secret('foo')
        then:
        config.getSecret() == ['foo']

        when:
        config.secret('bar')
        then:
        config.secret == ['foo', 'bar']
        config.getSecret() == ['foo', 'bar']
    }

    def 'should set process labels'() {
        when:
        def config = new ProcessConfig([:])
        then:
        config.getLabels() == []

        when:
        config.label('foo')
        then:
        config.getLabels() == ['foo']

        when:
        config.label('bar')
        then:
        config.getLabels() == ['foo','bar']
    }

    def 'should apply resource labels config' () {
        given:
        def config = new ProcessConfig(Mock(BaseScript))
        expect:
        config.getResourceLabels() == [:]

        when:
        config.resourceLabels([foo: 'one', bar: 'two'])
        then:
        config.getResourceLabels() == [foo: 'one', bar: 'two']

        when:
        config.resourceLabels([foo: 'new one', baz: 'three'])
        then:
        config.getResourceLabels() == [foo: 'new one', bar: 'two', baz: 'three']

    }

    def 'should check a valid label' () {

        expect:
        new ProcessConfig([:]).isValidLabel(lbl) == result

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
        ProcessConfig.matchesSelector(TARGET, SELECTOR) == EXPECTED

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
        def process = new ProcessConfig([:])
        process.applyConfigSelectorWithLabels(settings, ['short'])
        then:
        process.cpus == 1
        process.time == '1h'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithLabels(settings, ['long'])
        then:
        process.cpus == 32
        process.queue == 'cn-long'
        process.size() == 2

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithLabels(settings, ['foo'])
        then:
        process.cpus == 2
        process.disk == '100GB'
        process.queue == 'cn-long'
        process.size() == 3

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithLabels(settings, ['bar'])
        then:
        process.cpus == 32
        process.disk == '100GB'
        process.queue == 'cn-long'
        process.size() == 3

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithLabels(settings, ['gpu-1'])
        then:
        process.cpus == 4
        process.queue == 'cn-long'
        process.size() == 2

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
        def process = new ProcessConfig([:])
        process.applyConfigSelectorWithName(settings, 'xx')
        then:
        process.size() == 0

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithName(settings, 'alpha')
        then:
        process.cpus == 1
        process.time == '1h'
        process.size() == 2

        when:
        process =  new ProcessConfig([:])
        process.applyConfigSelectorWithName(settings, 'delta')
        then:
        process.cpus == 2
        process.disk == '100GB'
        process.size() == 2

        when:
        process =  new ProcessConfig([:])
        process.applyConfigSelectorWithName(settings, 'gamma')
        then:
        process.disk == '100GB'
        process.size() == 1

        when:
        process = new ProcessConfig([:])
        process.applyConfigSelectorWithName(settings, 'omega_x')
        then:
        process.cpus == 4
        process.size() == 1
    }


    def 'should apply config process defaults' () {

        when:
        def process = new ProcessConfig(Mock(BaseScript))

        // set process specific settings
        process.queue = 'cn-el6'
        process.memory = '10 GB'

        // apply config defaults
        process.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest'
        )

        then:
        process.queue == 'cn-el6'
        process.container == 'ubuntu:latest'
        process.memory == '10 GB'
        process.cacheable == true



        when:
        process = new ProcessConfig(Mock(BaseScript))
        // set process specific settings
        process.container = null
        // apply process defaults
        process.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        then:
        process.queue == 'def-queue'
        process.container == null
        process.maxRetries == 5



        when:
        process = new ProcessConfig(Mock(BaseScript))
        // set process specific settings
        process.maxRetries = 10
        // apply process defaults
        process.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        then:
        process.queue == 'def-queue'
        process.container == 'ubuntu:latest'
        process.maxRetries == 10
    }


    def 'should apply pod configs' () {

        when:
        def process =  new ProcessConfig([:])
        process.applyConfigDefaults( pod: [secret: 'foo', mountPath: '/there'] )
        then:
        process.pod == [
                [secret: 'foo', mountPath: '/there']
        ]

        when:
        process =  new ProcessConfig([:])
        process.applyConfigDefaults( pod: [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ] )

        then:
        process.pod == [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ]

    }

    def 'should clone config object' () {

        given:
        def config = new ProcessConfig(Mock(BaseScript))

        when:
        config.queue 'cn-el6'
        config.container 'ubuntu:latest'
        config.memory '10 GB'
        config._in_val('foo')
        config._in_file('sample.txt')
        config._out_file('result.txt')

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.getInputs().size() == 2
        config.getOutputs().size() == 1

        when:
        def copy = config.clone()
        copy.queue 'long'
        copy.container 'debian:wheezy'
        copy.memory '5 GB'
        copy._in_val('bar')
        copy._out_file('sample.bam')

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
        def process = new ProcessConfig(Mock(BaseScript))

        when:
        process.accelerator 5
        then:
        process.accelerator == [limit: 5]

        when:
        process.accelerator request: 1, limit: 5, type: 'nvida'
        then:
        process.accelerator == [request: 1, limit: 5, type: 'nvida']

        when:
        process.accelerator 5, type: 'nvida'
        then:
        process.accelerator == [limit: 5, type: 'nvida']

        when:
        process.accelerator 1, limit: 5
        then:
        process.accelerator == [request: 1, limit:5]

        when:
        process.accelerator 5, request: 1
        then:
        process.accelerator == [request: 1, limit:5]
    }

    def 'should apply disk config' () {

        given:
        def process = new ProcessConfig(Mock(BaseScript))

        when:
        process.disk '100 GB'
        then:
        process.disk == [request: '100 GB']

        when:
        process.disk '375 GB', type: 'local-ssd'
        then:
        process.disk == [request: '375 GB', type: 'local-ssd']

        when:
        process.disk request: '375 GB', type: 'local-ssd'
        then:
        process.disk == [request: '375 GB', type: 'local-ssd']
    }

    def 'should apply architecture config' () {

        given:
        def process = new ProcessConfig(Mock(BaseScript))

        when:
        process.arch 'linux/x86_64'
        then:
        process.arch == [name: 'linux/x86_64']

        when:
        process.arch 'linux/x86_64', target: 'zen3'
        then:
        process.arch == [name: 'linux/x86_64', target: 'zen3']

        when:
        process.arch name: 'linux/x86_64', target: 'zen3'
        then:
        process.arch == [name: 'linux/x86_64', target: 'zen3']
    }

    def 'should apply resourceLimits' () {
        given:
        def process = new ProcessConfig(Mock(BaseScript))

        when:
        process.resourceLimits time:'1h', memory: '2GB'
        then:
        process.resourceLimits == [time:'1h', memory: '2GB']
    }


    def 'should get default config path' () {
        given:
        ProviderConfig.env.remove('NXF_SCM_FILE')

        when:
        def path = ProviderConfig.getScmConfigPath()
        then:
        path.toString() == "${System.getProperty('user.home')}/.nextflow/scm"

    }

    def 'should get custom config path' () {
        given:
        def cfg = Files.createTempFile('test','config')
        ProviderConfig.env.NXF_SCM_FILE = cfg.toString()

        when:
        def path = ProviderConfig.getScmConfigPath()
        then:
        path.toString() == cfg.toString()

        cleanup:
        ProviderConfig.env.remove('NXF_SCM_FILE')
        cfg.delete()
    }

    def 'should not apply config on negative label' () {
        given:
        def settings = [
                'withLabel:foo': [ cpus: 2 ],
                'withLabel:!foo': [ cpus: 4 ],
                'withLabel:!nodisk_.*': [ disk: '100.GB']
        ]

        when:
        def p1 = new ProcessConfig([label: ['foo', 'other']])
        p1.applyConfig(settings, "processName", null, null)
        then:
        p1.cpus == 2
        p1.disk == '100.GB'

        when:
        def p2 = new ProcessConfig([label: ['foo', 'other', 'nodisk_label']])
        p2.applyConfig(settings, "processName", null, null)
        then:
        p2.cpus == 2
        !p2.disk

        when:
        def p3 = new ProcessConfig([label: ['other', 'nodisk_label']])
        p3.applyConfig(settings, "processName", null, null)
        then:
        p3.cpus == 4
        !p3.disk

    }

    def 'should throw exception for invalid error strategy' () {
        when:
        def process1 = new ProcessConfig(Mock(BaseScript))
        process1.errorStrategy 'abort'

        then:
        def e1 = thrown(IllegalArgumentException)
        e1.message == "Unknown error strategy 'abort' ― Available strategies are: terminate,finish,ignore,retry"

    }

    def 'should not throw exception for valid error strategy or closure' () {
        when:
        def process1 = new ProcessConfig(Mock(BaseScript))
        process1.errorStrategy 'retry'

        then:
        def e1 = noExceptionThrown()

        when:
        def process2 = new ProcessConfig(Mock(BaseScript))
        process2.errorStrategy 'terminate'

        then:
        def e2 = noExceptionThrown()

        when:
        def process3 = new ProcessConfig(Mock(BaseScript))
        process3.errorStrategy { task.exitStatus==14 ? 'retry' : 'terminate' }

        then:
        def e3 = noExceptionThrown()
    }
}
