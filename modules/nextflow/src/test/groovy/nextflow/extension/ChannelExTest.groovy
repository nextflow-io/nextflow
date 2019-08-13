package nextflow.extension

import nextflow.Channel
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.script.ChannelOut
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelExTest extends Specification {


    def 'should close the dataflow channel'() {

        when:
        def source = Channel.create()
        source << 10
        source << 20
        source << 30
        def result = source.close()
        then:
        result.is source
        result.val == 10
        result.val == 20
        result.val == 30
        result.val == Channel.STOP

        when:
        source = Channel.value().close()
        then:
        source.val == Channel.STOP

        when:
        source = Channel.value(1).close()
        then:
        source.val == 1

    }


    def 'should assign multiple channels in the current binding' () {
        given:
        def session = new Session()
        def ch1 = Channel.value('X')
        def ch2 = Channel.value('Y')

        when:
         new ChannelOut([ch1 ])
                 .set { alpha }
        then:
        session.binding.alpha.val == 'X'

        when:
        new ChannelOut([ch1, ch2 ])
                .set { alpha }
        then:
        thrown(ScriptRuntimeException)

    }


}
