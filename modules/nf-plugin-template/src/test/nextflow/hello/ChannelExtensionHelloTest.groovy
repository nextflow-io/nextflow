package nextflow.hello

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Session
import spock.lang.Specification


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class ChannelExtensionHelloTest extends Specification{

    def "should create a channel from hello"(){

        given:
        def session = Mock(Session)

        and:
        def helloExtension = new HelloExtension(); helloExtension.init(session)

        when:
        def result = helloExtension.reverse("Hi")

        then:
        result.val == 'iH'
        result.val == Channel.STOP
    }

    def "should consume a message from script"(){

        given:
        def session = Mock(Session)

        and:
        def helloExtension = new HelloExtension(); helloExtension.init(session)

        and:
        def ch = new DataflowQueue()
        ch.bind('Goodbye folks')
        ch.bind( Channel.STOP )

        when:
        def result = helloExtension.goodbye(ch)

        then:
        result.val == 'Goodbye folks'
        result.val == Channel.STOP
        helloExtension.goodbyeMessage == 'Goodbye folks'.toUpperCase()
    }
}
