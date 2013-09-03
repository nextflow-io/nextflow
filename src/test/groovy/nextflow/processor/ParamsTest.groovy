package nextflow.processor

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ParamsTest extends Specification {

    def testInParam()  {

        setup:
        def var = new DataflowVariable()
        var << 1

        when:
        def param = new FileInParam(name:'val', channel: var)
        then:
        param.name == 'val'
        param.channel.val == 1

        when:
        param = new StdInParam(channel: var)
        then:
        param.name == '-'
        param.channel.val == 1

    }


    def testOutParam() {

        when:
        def param = new FileOutParam(name:'val', channel: new DataflowVariable())
        then:
        param.name == 'val'
        param.channel instanceof DataflowVariable

        when:
        param = new StdOutParam(channel: new DataflowVariable())
        then:
        param.name == '-'
        param.channel instanceof DataflowVariable

    }

    def testCreateInParam() {

        /*
         * File input param def
         */
        when:
        def input = InParam.create( file: 'name.fa', from: 'hola' )
        then:
        input instanceof FileInParam
        input.name == 'name.fa'
        input.channel instanceof DataflowVariable
        input.channel.val == 'hola'


        /*
         * 'stdin' input param
         */
        when:
        input = InParam.create( file: '-', from: 'hola' )
        then:
        input instanceof StdInParam
        input.name == '-'
        input.channel instanceof DataflowVariable
        input.channel.val == 'hola'

        /*
         * 'value' input param
         */
        when:
        input = InParam.create( val: 'ref_name', from: [1,2,3] )
        then:
        input instanceof ValueInParam
        input.name == 'ref_name'
        input.channel instanceof DataflowQueue

        /*
         * 'environment' input param
         */
        when:
        input = InParam.create( env: 'var_name', from: 10 )
        then:
        input instanceof EnvInParam
        input.name == 'var_name'
        input.channel instanceof DataflowVariable

        // missing type
        when:
        InParam.create( from: 10 )
        then:
        thrown(IllegalArgumentException)

    }


    def testCreateOutParam() {

        setup:
        def channel3 = new DataflowVariable()
        channel3 << 10

        def script = new Script() {
            @Override
            Object run() { return null }
        }
        script.binding.setVariable('channel3', channel3)

        /*
         * test that:
         * - create a 'FileOutParam'
         * - instantiate a DataflowQueue as channel
         * - the channel is bound to the script variables context
         */
        when:
        def out = OutParam.create( file:'simple.fa', into: 'channel1', script )
        then:
        out instanceof FileOutParam
        out.name == 'simple.fa'
        out.channel instanceof DataflowQueue
        script.binding.getVariable('channel1') instanceof DataflowQueue

        /*
         * test that:
         * - create a 'StdOutParam' since the special name '-' is used
         * - instantiate a DataflowQueue as channel
         * - the channel is *NOT* bound to the script variables context
         */
        when:
        out = OutParam.create( file:'-', into: 'channel2', script )
        then:
        out instanceof StdOutParam
        out.name == '-'
        out.channel instanceof DataflowQueue
        script.binding.hasVariable('channel2')


        /*
         * test that:
         * - use the channel instance defined in the script bindings
         */
        when:
        out = OutParam.create( file:'file.txt', into: 'channel3', script )
        then:
        out.name == 'file.txt'
        out.channel instanceof DataflowVariable
        script.binding.getVariable('channel3').is( channel3 )

        /*
         *
         */
        when:
        def out1 = OutParam.create( file:'file.txt', into: 'channel4', autoClose: false, joint: true, script )
        def out2 = OutParam.create( file:'file.txt', into: 'channel4', autoClose: true, joint: false, script )

        then:
        !(out1 as FileOutParam).autoClose
        (out1 as FileOutParam).joint

        (out2 as FileOutParam).autoClose
        !(out2 as FileOutParam).joint



    }

}
