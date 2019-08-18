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
class CHTest extends Specification {

    def 'should create dataflow variable or queue' () {

        expect:
        CH.create() instanceof DataflowQueue
        CH.create(false) instanceof DataflowQueue
        CH.create(true) instanceof DataflowVariable

        CH.createBy(new DataflowVariable()) instanceof DataflowVariable
        CH.createBy(new DataflowQueue()) instanceof DataflowQueue


        when:
        NextflowMeta.instance.enableDsl2()
        then:
        CH.create() instanceof DataflowBroadcast
        CH.create(false) instanceof DataflowBroadcast
        CH.create(true) instanceof DataflowVariable

        CH.createBy(new DataflowVariable()) instanceof DataflowVariable
        CH.createBy(new DataflowQueue()) instanceof DataflowBroadcast

        cleanup:
        NextflowMeta.instance.disableDsl2()

    }


    def 'should check queue channel' () {
        expect:
        CH.isChannelQueue(new DataflowQueue())
        CH.isChannelQueue(new DataflowBroadcast().createReadChannel())
        !CH.isChannelQueue(new DataflowVariable())
        !CH.isChannelQueue('hello')
    }

    def 'should validate allScalar method' () {

        expect:
        CH.allScalar([1])
        CH.allScalar([1, 2, 3])
        !CH.allScalar([new DataflowVariable(), new DataflowQueue()])
        !CH.allScalar([new DataflowQueue(), new DataflowQueue()])
    }

}
