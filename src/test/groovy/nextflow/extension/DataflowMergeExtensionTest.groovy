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
class DataflowMergeExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }

    def 'should merge with open array'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);
        def delta = Channel.from(7,8,1);

        def result = alpha.merge( beta, delta ) { a,b,c -> [a,b,c] }

        then:
        result instanceof DataflowQueue
        result.val == [1,2,7]
        result.val == [3,4,8]
        result.val == [5,6,1]
        result.val == Channel.STOP
    }

    def 'should merge with list'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);
        def delta = Channel.from(7,8,1);

        def result = alpha.merge( [beta, delta] ) { a,b,c -> [c,b,a] }

        then:
        result instanceof DataflowQueue
        result.val == [7,2,1]
        result.val == [8,4,3]
        result.val == [1,6,5]
        result.val == Channel.STOP
    }

    def 'should merge with queue'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);

        def result = alpha.merge(beta) { a,b -> [a, b+1] }

        then:
        result instanceof DataflowQueue
        result.val == [1,3]
        result.val == [3,5]
        result.val == [5,7]
        result.val == Channel.STOP
    }

    def 'should merge with variables'() {

        when:
        def alpha = Channel.value('Hello');
        def beta = Channel.value('World')

        def result = alpha.merge(beta) { a,b -> [a, b] }

        then:
        result instanceof DataflowVariable
        result.val == ['Hello','World']
        result.val == ['Hello','World']
    }

}
