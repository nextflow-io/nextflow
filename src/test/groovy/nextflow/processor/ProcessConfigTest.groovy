/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.processor
import static nextflow.util.CacheHelper.HashMode

import nextflow.exception.IllegalDirectiveException
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
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
        config.validExitStatus == [0]
        config.maxRetries == 0
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
        def config = new ProcessConfig(script).enterCaptureMode(true)
        def x = config.hola

        then:
        thrown(MissingPropertyException)
    }


    def 'should check property existence' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        expect:
        config.containsKey('echo')
        config.containsKey('shell')
        config.containsKey('validExitStatus')
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
        config._in_val('x') .from(1)
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
        config._out_file(new TokenVar('file1')).into('ch1')
        config._out_file(new TokenVar('file2')).into('ch2')
        config._out_file(new TokenVar('file3')).into('ch3')

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
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config.publishDir '/data'
        then:
        config.get('publishDir') == '/data'

        when:
        config.publishDir '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir') == [[mode: 'link', pattern: '*.bam'], '/data']

        when:
        config.publishDir path: '/data', mode: 'link', pattern: '*.bam'
        then:
        config.get('publishDir') == [path: '/data', mode: 'link', pattern: '*.bam']
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

}
