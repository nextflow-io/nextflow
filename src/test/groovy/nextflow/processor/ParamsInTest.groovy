package nextflow.processor
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Nextflow
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsInTest extends Specification {



    def testCreateInParam() {

        setup:
        def binding = new Binding()
        def script = Mock(Script)
        script.getBinding() >> binding

        /*
         * File input parameter
         * - when passing with *using* a ref channel, it will ne the one used to get the value
         */
        when:
        def var = Nextflow.val(1)
        def param = new FileInParam(script, 'x').using(var)
        then:
        param.name == 'x'
        param.channel.val == 1


        /*
         * File input param def
         * - when passing with *using* a generic value NOT a channel, a new channel instance
         *   it is created and associated to the parameter
         */
        when:
        def input = new FileInParam(script, 'name.fa' ) .using('hola')
        then:
        input instanceof FileInParam
        input.name == 'name.fa'
        input.channel instanceof DataflowVariable
        input.channel.val == 'hola'


        /*
         * 'stdin' input param
         * - when passing with *using* a generic value NOT a channel, a new channel instance
         *   it is created and associated to the parameter
         */
        when:
        input = new StdInParam(script).using('hola')
        then:
        input instanceof StdInParam
        input.name == '-'
        input.channel instanceof DataflowVariable
        input.channel.val == 'hola'

        /*
         * 'value' input param
         */
        when:
        input = new ValueInParam(script, 'ref_name' ).using( [1,2,3] )
        then:
        input instanceof ValueInParam
        input.name == 'ref_name'
        input.channel instanceof DataflowQueue
        input.channel.val == 1
        input.channel.val == 2
        input.channel.val == 3

        /*
         * 'environment' input param
         */
        when:
        input = new EnvInParam(script, 'var_name' ).using(10)
        then:
        input instanceof EnvInParam
        input.name == 'var_name'
        input.channel instanceof DataflowVariable


        /*
         * test access to script context variable
         *
         * - when the channel is NOT specified by *using* method
         *   it tries to access the script variable having the same name specified
         *
         */
        when:
        binding.setVariable('a_script_value', Nextflow.val(3) )
        input = new ValueInParam(script, 'a_script_value' )
        then:
        input.name == 'a_script_value'
        input.channel == script.getBinding().getVariable('a_script_value' )


        /*
         * test access to script context variable
         *
         * - like before access the value in the script context
         * - since it is not a dataflow channel, wrap the value by a new channel instance
         */
        when:
        binding.setVariable('a_script_x', 4 )
        input = new ValueInParam(script, 'a_script_x' )
        then:
        input.name == 'a_script_x'
        input.channel instanceof DataflowVariable
        input.channel.val == 4


        /*
         * - No channel has been specified by *using* method
         * - no variable exist in the script context with variable name
         * - it results to an exception
         */
        when:
        input = new ValueInParam(script, 'a_script_z' )
        def x = input.getChannel()
        then:
        thrown(IllegalStateException)

    }


}
