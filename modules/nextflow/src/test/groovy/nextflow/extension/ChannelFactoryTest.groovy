package nextflow.extension

import spock.lang.Specification

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.NextflowMeta
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelFactoryTest extends Specification {

    def 'should create dataflow variable or queue' () {

        expect:
        ChannelFactory.create() instanceof DataflowQueue
        ChannelFactory.create(false) instanceof DataflowQueue
        ChannelFactory.create(true) instanceof DataflowVariable

        ChannelFactory.createBy(new DataflowVariable()) instanceof DataflowVariable
        ChannelFactory.createBy(new DataflowQueue()) instanceof DataflowQueue


        when:
        NextflowMeta.instance.enableDsl2()
        then:
        ChannelFactory.create() instanceof DataflowBroadcast
        ChannelFactory.create(false) instanceof DataflowBroadcast
        ChannelFactory.create(true) instanceof DataflowVariable

        ChannelFactory.createBy(new DataflowVariable()) instanceof DataflowVariable
        ChannelFactory.createBy(new DataflowQueue()) instanceof DataflowBroadcast

        cleanup:
        NextflowMeta.instance.disableDsl2()

    }


    def 'should check queue channel' () {
        expect:
        ChannelFactory.isChannelQueue(new DataflowQueue())
        ChannelFactory.isChannelQueue(new DataflowBroadcast().createReadChannel())
        !ChannelFactory.isChannelQueue(new DataflowVariable())
        !ChannelFactory.isChannelQueue('hello')
    }

    def 'should validate allScalar method' () {

        expect:
        ChannelFactory.allScalar([1])
        ChannelFactory.allScalar([1, 2, 3])
        !ChannelFactory.allScalar([new DataflowVariable(), new DataflowQueue()])
        !ChannelFactory.allScalar([new DataflowQueue(), new DataflowQueue()])
    }

}
