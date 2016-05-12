package nextflow.extension

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowTapExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }


    def 'should `tap` item to a new channel' () {

        when:
        def result = Channel.from( 4,7,9 ) .tap { first }.map { it+1 }
        then:
        session.binding.first.val == 4
        session.binding.first.val == 7
        session.binding.first.val == 9
        session.binding.first.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` target channel' () {

        when:
        def target = Channel.create()
        def result = Channel.from( 8,2,5 ) .tap(target).map { it+1 }
        then:
        result instanceof DataflowQueue
        target instanceof DataflowQueue

        target.val == 8
        target.val == 2
        target.val == 5
        target.val == Channel.STOP

        result.val == 9
        result.val == 3
        result.val == 6
        result.val == Channel.STOP

        !session.dag.isEmpty()

    }

    def 'should `tap` dataflow value' () {

        when:
        def target = Channel.value()
        def result = Channel.value(7) .tap(target).map { it+1 }
        then:
        result instanceof DataflowVariable
        target instanceof DataflowVariable

        target.val == 7
        target.val == 7
        result.val == 8
        result.val == 8

        !session.dag.isEmpty()

    }

    def 'should `tap` dataflow value and target as queue' () {

        when:
        new DataflowVariable() .tap( new DataflowQueue() )
        then:
        thrown(IllegalArgumentException)
        session.dag.isEmpty()

    }

}
