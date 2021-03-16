package nextflow.script

import nextflow.Channel
import nextflow.script.params.BaseOutParam
import nextflow.script.params.OutputsList
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestflowDslTest extends Specification {

    def 'should create testflow output'() {
        given:
        def queue1 = Channel.of(1,2,3)
        def queue2 = Channel.of('A', 'B', 'C')
        def output = new ChannelOut([queue1, queue2])

        when:
        def flow = new TestflowDsl(output)
        then:
        flow.emissionsCount() == 3
        and:
        flow.emissionNext {
            assert out[0] == 1
            assert out[1] == 'A'
        }
        and:
        flow.emissionNext {
            assert out[0] == 2
            assert out[1] == 'B'
        }
        and:
        flow.emissionNext {
            assert out[0] == 3
            assert out[1] == 'C'
        }
    }

    def 'should create testflow output 2'() {
        given:
        def queue1 = Channel.of(1,2)
        def queue2 = Channel.of('A', 'B', 'C', 'D')
        def output = new ChannelOut([queue1, queue2])

        when:
        def flow = new TestflowDsl(output)
        then:
        flow.emissionsCount() == 2
        and:
        flow.emissionNext {
            assert out[0] == 1
            assert out[1] == 'A'
        }
        and:
        flow.emissionNext {
            assert out[0] == 2
            assert out[1] == 'B'
        }
    }

    def 'should create testflow output 3'() {
        given:
        def queue1 = Channel.value(10)
        def queue2 = Channel.of('A', 'B', 'C')
        def output = new ChannelOut([queue1, queue2])

        when:
        def flow = new TestflowDsl(output)
        then:
        flow.emissionsCount() == 3
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[1] == 'A'
        }
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[1] == 'B'
        }
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[1] == 'C'
        }
    }

    def 'should create testflow with named params' () {
        given:
        def queue1 = Channel.value(10)
        def queue2 = Channel.of('A', 'B', 'C')
        and:
        def param1 = Mock(BaseOutParam) {
            getChannelEmitName() >> 'foo'
            getOutChannel() >> queue1
        }
        and:
        def param2 = Mock(BaseOutParam) {
            getChannelEmitName() >> 'bar'
            getOutChannel() >> queue2
        }
        and:
        def params = new OutputsList()
        params << param1
        params << param2

        when:
        def output = new ChannelOut(params)
        def flow = new TestflowDsl(output)

        then:
        flow.emissionsCount() == 3
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[0] == out.foo
            assert out[1] == 'A'
            assert out[1] == out.bar
        }
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[1] == 'B'
            assert out[0] == out.foo
            assert out[1] == out.bar
        }
        and:
        flow.emissionNext{
            assert out[0] == 10
            assert out[1] == 'C'
            assert out[0] == out.foo
            assert out[1] == out.bar
        }
    }

    def 'should create testflow with all value' () {
        given:
        def queue1 = Channel.value(10)
        def queue2 = Channel.value(20)
        and:
        def param1 = Mock(BaseOutParam) {
            getChannelEmitName() >> 'foo'
            getOutChannel() >> queue1
        }
        and:
        def param2 = Mock(BaseOutParam) {
            getChannelEmitName() >> 'bar'
            getOutChannel() >> queue2
        }
        and:
        def params = new OutputsList()
        params << param1
        params << param2

        when:
        def output = new ChannelOut(params)
        def flow = new TestflowDsl(output)

        then:
        flow.emissionsCount() == 1
        and:
        flow.emissionNext {
            assert out[0] == 10
            assert out[0] == out.foo
            assert out[1] == 20
            assert out[1] == out.bar
        }

    }
}
