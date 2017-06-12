package nextflow.extension

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowHelperTest extends Specification {


    def setupSpec() {
        new Session()
    }

    def testHandlerNames() {

        when:
        DataflowHelper.checkSubscribeHandlers( [:] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}] )
        then:
        true

        when:
        DataflowHelper.checkSubscribeHandlers( [ onNext:{}, xxx:{}] )
        then:
        thrown(IllegalArgumentException)

        when:
        DataflowHelper.checkSubscribeHandlers( [ xxx:{}] )
        then:
        thrown(IllegalArgumentException)
    }

}
