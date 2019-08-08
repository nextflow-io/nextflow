package nextflow.extension

import nextflow.Session
import nextflow.script.ChannelArrayList
import spock.lang.Specification

import nextflow.Channel

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
        def ch3 = Channel.value('Z')
        def output = new ChannelArrayList([ch1, ch2, ch3])

        when:
        output.set { alpha; bravo; delta }

        then:
        session.binding.alpha.val == 'X'
        session.binding.bravo.val == 'Y'
        session.binding.delta.val == 'Z'

        // should throw an exception because
        // defines more channel variables
        // then existing ones
        when:
        output.set { X; Y; W; Z }
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Operation `set` expects 4 channels but only 3 are provided"

        when:
        output.set { ALPHA; ALPHA; BRAVO }
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'Duplicate channel definition: ALPHA'
    }


}
