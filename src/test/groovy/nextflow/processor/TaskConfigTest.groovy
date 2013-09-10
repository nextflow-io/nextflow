package nextflow.processor

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


    def 'test defaults' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

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
        config.containsKey('validExitCodes')
        config.containsKey('inputs')
        config.containsKey('outputs')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null


    }


    def 'test input' () {

        setup:
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        when:
        config.__in_file('filename.fa') .using(new DataflowVariable<>())
        config.__in_val('x') .using(1)
        config.stdin().using(new DataflowVariable<>())

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
        def script = Mock(BaseScript)
        def config = new TaskConfig(script)

        when:
        config.stdout()
        config.__out_file('file1.fa').using('ch1')
        config.__out_file('file2.fa').using('ch2')
        config.__out_file('file3.fa').using('ch3')

        then:
        config.outputs.size() == 4
        config.outputs.names == ['-', 'file1.fa', 'file2.fa', 'file3.fa']
        config.outputs.ofType(StdOutParam).size() == 1

        config.outputs[0] instanceof StdOutParam
        config.outputs[1].name == 'file1.fa'
        config.outputs[2].name == 'file2.fa'
        config.outputs[3].name == 'file3.fa'


    }


}
