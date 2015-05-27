/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.script.ValueSharedParam
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
    }

    def 'should set properties' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        // setting property using method without brackets
        when:
        config.hola 'val 1'
        then:
        config.hola == 'val 1'

        // setting list values
        when:
        config.hola 1,2,3
        then:
        config.hola == [1,2,3]

        // setting named parameters attribute
        when:
        config.hola field1:'val1', field2: 'val2'
        then:
        config.hola == [field1:'val1', field2: 'val2']

        // generic value assigned like a 'plain' property
        when:
        config.hola = 99
        then:
        config.hola == 99

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
        config._out_file('file1.fa').into('ch1')
        config._out_file('file2.fa').into('ch2')
        config._out_file('file3.fa').into('ch3')

        then:
        config.outputs.size() == 4
        config.outputs.names == ['-', 'file1.fa', 'file2.fa', 'file3.fa']
        config.outputs.ofType(StdOutParam).size() == 1

        config.outputs[0] instanceof StdOutParam
        config.outputs[1].name == 'file1.fa'
        config.outputs[2].name == 'file2.fa'
        config.outputs[3].name == 'file3.fa'

    }

    /*
     *  shared: val x (seed x)
     *  shared: val x seed y
     *  shared: val x seed y using z
     */
    def testSharedValue() {

        setup:
        def binding = new Binding()
        def script = Mock(BaseScript)
        script.getBinding() >> { binding }

        when:
        def config = new ProcessConfig(script)
        def val = config._share_val( new TokenVar('xxx'))
        then:
        val instanceof ValueSharedParam
        val.name == 'xxx'
        val.inChannel.val == null
        val.outChannel == null

        when:
        binding.setVariable('yyy', 'Hola')
        config = new ProcessConfig(script)
        val = config._share_val(new TokenVar('yyy'))
        then:
        val instanceof ValueSharedParam
        val.name == 'yyy'
        val.inChannel.val == 'Hola'
        val.outChannel == null

        // specifying a value with the 'using' method
        // that value is bound to the input channel
        when:
        config = new ProcessConfig(script)
        val = config._share_val('yyy') .from('Beta')
        then:
        val instanceof ValueSharedParam
        val.name == 'yyy'
        val.inChannel.val == 'Beta'
        val.outChannel == null

        // specifying a 'closure' with the 'using' method
        // that value is bound to the input channel
        when:
        config = new ProcessConfig(script)
        val = config._share_val('yyy') .from({ 99 })
        then:
        val instanceof ValueSharedParam
        val.name == 'yyy'
        val.inChannel.val == 99
        val.outChannel == null


        // specifying a 'channel' it is reused
        // that value is bound to the input channel
        when:
        def channel = new DataflowVariable()
        channel << 123

        config = new ProcessConfig(script)
        val = config._share_val('zzz') .from(channel)
        then:
        val instanceof ValueSharedParam
        val.name == 'zzz'
        val.inChannel.getVal() == 123
        val.outChannel == null

        // when a channel name is specified with the method 'into'
        // a DataflowVariable is created in the script context
        when:
        config = new ProcessConfig(script)
        val = config._share_val(new TokenVar('x1')) .into( new TokenVar('x2') )
        then:
        val instanceof ValueSharedParam
        val.name == 'x1'
        val.inChannel.getVal() == null
        val.outChannel instanceof DataflowVariable
        binding.getVariable('x2') == val.outChannel

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

}
