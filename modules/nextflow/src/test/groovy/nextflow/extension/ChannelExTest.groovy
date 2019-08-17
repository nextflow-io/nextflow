package nextflow.extension

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


}
