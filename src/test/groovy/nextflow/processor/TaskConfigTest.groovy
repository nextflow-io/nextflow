package nextflow.processor

import groovyx.gpars.dataflow.DataflowVariable
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import nextflow.script.InputsList
import nextflow.script.OutputsList
import nextflow.script.ScriptVar
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.ValueInParam
import nextflow.script.ValueSharedParam
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskConfigTest extends Specification {


    def 'test defaults' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        expect:
        config.shell ==  ['/bin/bash','-ue']
        config.cacheable
        config.validExitStatus == [0]
        config.errorStrategy == ErrorStrategy.TERMINATE
        config.inputs instanceof InputsList
        config.outputs instanceof OutputsList

    }


    def 'test merge' ()  {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        expect:
        !config.merge

        when:
        config.setProperty('merge',true)
        then:
        config.merge



    }

    def 'test setting properties' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

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

        // maxDuration property
        when:
        config.maxDuration '1h'
        then:
        config.maxDuration == new Duration('1h')
        config.maxDuration as Duration == new Duration('1h')

        // maxMemory property
        when:
        config.maxMemory '2GB'
        then:
        config.maxMemory == new MemoryUnit('2GB')

        // generic value assigned like a 'plain' property
        when:
        config.hola = 99
        then:
        config.hola == 99

    }

    def testParseProperties() {

        when:
        def config = new TaskConfig( maxDuration:'1h' )
        then:
        config.maxDuration as Duration == Duration.of('1h')
    }


    def 'test NO missingPropertyException' () {

        when:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)
        def x = config.hola

        then:
        x == null
        noExceptionThrown()

    }

    def 'test MissingPropertyException' () {
        when:
        def script = Mock(BaseScript)
        def config = new TaskConfigWrapper(new TaskConfig(script))
        def x = config.hola

        then:
        thrown(MissingPropertyException)
    }


    def 'test check property existence' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        expect:
        config.containsKey('echo')
        config.containsKey('shell')
        config.containsKey('validExitStatus')
        config.containsKey('inputs')
        config.containsKey('outputs')
        config.containsKey('undef')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null

    }

    def 'test undef' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        expect:
        config.undef == false

        when:
        config.undef(true)
        then:
        config.undef == true

    }


    def 'test input' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

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

    def 'test outputs' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

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
        def config = new TaskConfig(script)
        def val = config._share_val( new ScriptVar('xxx'))
        then:
        val instanceof ValueSharedParam
        val.name == 'xxx'
        val.inChannel.val == null
        val.outChannel == null

        when:
        binding.setVariable('yyy', 'Hola')
        config = new TaskConfig(script)
        val = config._share_val(new ScriptVar('yyy'))
        then:
        val instanceof ValueSharedParam
        val.name == 'yyy'
        val.inChannel.val == 'Hola'
        val.outChannel == null

        // specifying a value with the 'using' method
        // that value is bound to the input channel
        when:
        config = new TaskConfig(script)
        val = config._share_val('yyy') .from('Beta')
        then:
        val instanceof ValueSharedParam
        val.name == 'yyy'
        val.inChannel.val == 'Beta'
        val.outChannel == null

        // specifying a 'closure' with the 'using' method
        // that value is bound to the input channel
        when:
        config = new TaskConfig(script)
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

        config = new TaskConfig(script)
        val = config._share_val('zzz') .from(channel)
        then:
        val instanceof ValueSharedParam
        val.name == 'zzz'
        val.inChannel.getVal() == 123
        val.outChannel == null

        // when a channel name is specified with the method 'into'
        // a DataflowVariable is created in the script context
        when:
        config = new TaskConfig(script)
        val = config._share_val(new ScriptVar('x1')) .into( new ScriptVar('x2') )
        then:
        val instanceof ValueSharedParam
        val.name == 'x1'
        val.inChannel.getVal() == null
        val.outChannel instanceof DataflowVariable
        binding.getVariable('x2') == val.outChannel

    }


}
