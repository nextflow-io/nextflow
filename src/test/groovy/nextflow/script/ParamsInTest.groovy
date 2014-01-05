package nextflow.script
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
        FileInParam param = new FileInParam(script, var)._as('x')

        then:
        param.name == null
        param.filePattern == 'x'
        param.inTarget == var
        param.inChannel.val == 1


        /*
         * File input param def
         * - when passing with *using* a generic value NOT a channel, a new channel instance
         *   it is created and associated to the parameter
         */
        when:
        def input = new FileInParam(script, 'hola')._as('name.fa')

        then:
        input instanceof FileInParam
        input.name == null
        input.filePattern == 'name.fa'
        input.inChannel instanceof DataflowVariable
        input.inChannel.val == 'hola'


        when:
        binding.setVariable('channel_x', Nextflow.val(1))
        input = new FileInParam(script, new ScriptVar('channel_x'))
        def x = input.getInChannel()
        then:
        input.name == 'channel_x'
        // TODO !!! filePattern returning '*'
        // ????????
        //input.filePattern == '*'
        input.inChannel == x
        input.inChannel.val == 1


        /*
         * 'stdin' input param
         * - when passing with *using* a generic value NOT a channel, a new channel instance
         *   it is created and associated to the parameter
         */
        when:
        input = new StdInParam(script,'hola')
        then:
        input instanceof StdInParam
        input.name == '-'
        input.inChannel instanceof DataflowVariable
        input.inChannel.val == 'hola'

        /*
         * 'value' input param
         */
        when:
        input = new ValueInParam(script, [1,2,3] )._as( 'ref_name' )
        then:
        input instanceof ValueInParam
        input.name == 'ref_name'
        input.inChannel instanceof DataflowQueue
        input.inChannel.val == 1
        input.inChannel.val == 2
        input.inChannel.val == 3

        /*
         * 'environment' input param
         */
        when:
        input = new EnvInParam(script, 10 )._as('var_name')
        then:
        input instanceof EnvInParam
        input.name == 'var_name'
        input.inChannel instanceof DataflowVariable


        /*
         * test access to script context variable
         *
         * - when the channel is NOT specified by *using* method
         *   it tries to access the script variable having the same name specified
         *
         */
        when:
        binding.setVariable('a_script_value', Nextflow.val(3) )
        input = new ValueInParam(script, new ScriptVar('a_script_value') )
        then:
        input.name == 'a_script_value'
        input.inChannel == script.getBinding().getVariable('a_script_value' )


        /*
         * test access to script context variable
         *
         * - like before access the value in the script context
         * - since it is not a dataflow channel, wrap the value by a new channel instance
         */
        when:
        binding.setVariable('a_script_x', 4 )
        input = new ValueInParam(script, new ScriptVar('a_script_x') )
        then:
        input.name == 'a_script_x'
        input.inChannel instanceof DataflowVariable
        input.inChannel.val == 4


        /*
         * - No channel has been specified by *using* method
         * - no variable exist in the script context with variable name
         * - it results to an exception
         */
        when:
        new ValueInParam(script, new ScriptVar('a_script_z' ))
        then:
        thrown(MissingPropertyException)

    }


}
