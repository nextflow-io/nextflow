package nextflow.processor

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.script.BaseScript
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskConfigTest extends Specification {

    static class DummyScript extends BaseScript {
        @Override
        Object run() {
            return null
        }
    }


    def 'test defaults' () {

        setup:
        def config = new TaskConfig(new DummyScript())

        expect:
        config.shell ==  ['/bin/bash','-ue']
        config.cacheable
        config.validExitCodes == [0]
        config.errorStrategy == ErrorStrategy.TERMINATE
        config.inputs instanceof InputsList
        config.outputs instanceof OutputsList

    }

    def 'test setting properties' () {

        setup:
        def config = new TaskConfig(new DummyScript())

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


    def 'test NO missingPropertyException' () {

        when:
        def config = new TaskConfig(new DummyScript())
        def x = config.hola

        then:
        x == null
        noExceptionThrown()

    }

    def 'test MissingPropertyException' () {
        when:
        def config = new TaskConfigWrapper(new TaskConfig(new DummyScript()))
        def x = config.hola

        then:
        thrown(MissingPropertyException)
    }


    def 'test check property existence' () {

        setup:
        def config = new TaskConfig(new DummyScript())

        expect:
        config.containsKey('echo')
        config.containsKey('shell')
        config.containsKey('validExitCodes')
        config.containsKey('inputs')
        config.containsKey('outputs')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null


    }


    def 'test input' () {

        setup:
        def config = new TaskConfig(new DummyScript())

        when:
        config.input file: 'filename.fa', from: new DataflowVariable<>()
        config.input val: 'x', from: 1
        config.input file:'-', from: new DataflowVariable<>()

        then:
        config.getInputs().size() == 3

        config.inputs.get(0) instanceof FileInParam
        config.inputs.get(0).name == 'filename.fa'

        config.inputs.get(1) instanceof ValueInParam
        config.inputs.get(1).name == 'x'

        config.inputs.get(2).name == '-'
        config.inputs.get(2) instanceof StdInParam

        config.inputs.names == [ 'filename.fa', 'x', '-' ]
        config.inputs.ofType( FileInParam ) == [ config.getInputs().get(0) ]

    }

    def 'test outputs' () {

        setup:
        def ch1 = new DataflowVariable()
        def ch3 = new DataflowVariable()
        ch1 << 1

        def script = new DummyScript()
        def config = new TaskConfig(script)
        script.setProperty('ch1', ch1)


        when:
        config.stdout( new DataflowVariable() )
        config.output file: 'file1.fa', into: 'ch1'
        config.output file: 'file2.fa', into: 'ch2'
        config.output file: 'file3.fa', into: ch3

        then:
        config.outputs.size() == 4
        config.outputs.names == ['-', 'file1.fa', 'file2.fa', 'file3.fa']
        config.outputs.ofType(StdOutParam).size() == 1

        config.outputs[0] instanceof StdOutParam
        config.outputs[1].name == 'file1.fa'
        config.outputs[2].name == 'file2.fa'
        config.outputs[3].name == 'file3.fa'

        // 'ch1' is the same channel defined in the script binding
        script.getBinding().getVariable('ch1') == ch1
        script.getBinding().getVariable('ch1').val == 1

        // 'ch2' is created automatically and added to the script context as a 'DataflowQueue'
        script.getBinding().getVariable('ch2') instanceof DataflowQueue

        // 'ch3' is not in the script context because it has specified by 'ref' (not by name)
        !script.getBinding().hasVariable('ch3')

    }


}
