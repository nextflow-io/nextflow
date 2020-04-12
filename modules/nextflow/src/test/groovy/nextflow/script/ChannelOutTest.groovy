package nextflow.script

import spock.lang.Specification

import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.exception.DuplicateChannelNameException
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelOutTest extends Specification {
    

    def 'should get out by name' () {
        given:
        def ch1 = Channel.value('a')
        def ch2 = Channel.value('b')

        when:
        def out = new ChannelOut([foo:ch1, bar:ch2])
        then:
        out[0].val == 'a'
        out[1].val == 'b'
        and:
        out.foo.val == 'a'
        out.bar.val == 'b'
        and:
        out.getNames() == ['foo','bar'] as Set
    }

    def 'should validate output spread' () {
        given:
        def out1 = new ChannelOut( ['a', 'b'] )
        def out2 = new ChannelOut( ['x', 'y', 'z'])

        expect:
        ChannelOut.spread([1, 2, 3]) \
            == [1, 2, 3]

        ChannelOut.spread([out1]) \
            == ['a', 'b']

        ChannelOut.spread([out1, 'p', 'q', out2] ) \
            == ['a', 'b', 'p', 'q', 'x', 'y', 'z']
    }

    def 'should check spread as array' () {
        given:
        def out1 = new ChannelOut( ['a', 'b'] )
        def out2 = new ChannelOut( ['x', 'y', 'z'])
        and:
        def ARGS1 = [1, 2, 3] as Object[]
        def ARGS2 = [out1] as Object[]
        def ARGS3 = [out1, 'p', 'q', out2] as Object[]

        expect:
        ChannelOut.spreadToArray(ARGS1) == ARGS1
        ChannelOut.spreadToArray(ARGS1).is(ARGS1)
        and:
        ChannelOut.spreadToArray(ARGS2) == ['a', 'b'] as Object[]
        and:
        ChannelOut.spreadToArray(ARGS3) == ['a', 'b', 'p', 'q', 'x', 'y', 'z'] as Object[]
    }

    def 'should create with outputs list' () {
        given:
        def ch1 = Mock(DataflowWriteChannel)
        def ch2 = Mock(DataflowWriteChannel)
        def ch3 = Mock(DataflowWriteChannel)
        def ch4 = Mock(DataflowWriteChannel)

        def p1 = Mock(OutParam) { getOutChannel() >> ch1; getChannelEmitName() >> 'foo'; }
        def p2 = Mock(OutParam) { getOutChannel() >> ch2; getChannelEmitName() >> 'bar'; }
        def p3 = Mock(OutParam) { getOutChannel() >> ch3; getChannelEmitName() >> null }
        def p4 = Mock(OutParam) { getOutChannel() >> ch3; getChannelEmitName() >> 'foo' }

        def list = new OutputsList()
        list.add(p1)
        list.add(p2)
        list.add(p3)

        when:
        def out = new ChannelOut(list)
        then:
        out.size() == 3
        and:
        out.getNames() == ['foo', 'bar'] as Set
        out[0] == ch1
        out[1] == ch2
        out[2] == ch3
        out.foo == ch1
        out.bar == ch2

        when:
        list.add(p4)
        new ChannelOut(list)
        then:
        def e = thrown(DuplicateChannelNameException)
        e.message == 'Output channel name `foo` is used more than one time'

    }
}
