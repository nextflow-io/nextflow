package nextflow.extension

import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowHelperTest extends Specification {


    def setupSpec() {
        new Session()
    }

    def 'should subscribe handlers'() {

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

    @Unroll
    def 'should split entry' () {
        when:
        def pair = DataflowHelper.split(pivot, entry)
        then:
        pair.keys == keys
        pair.values == values

        where:
        pivot           | entry                         | keys          | values
        [0]             | ['A','B','C','D','E','F']     | ['A']         | ['B','C','D','E','F']
        [0,1]           | ['A','B','C','D','E','F']     | ['A','B']     | ['C','D','E','F']
        [0,2]           | ['A','B','C','D','E','F']     | ['A','C']     | ['B','D','E','F']
        [0,1,4]         | ['A','B','C','D','E','F']     | ['A','B','E'] | ['C','D','F']
        [0]             | 'A'                           | ['A']         | []
    }
}
