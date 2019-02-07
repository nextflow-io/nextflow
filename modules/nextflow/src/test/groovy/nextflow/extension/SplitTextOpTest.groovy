package nextflow.extension

import spock.lang.Specification

import nextflow.Channel

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitTextOpTest extends Specification {

    def 'should split text' () {

        when:
        def result = Channel.from('foo\nbar').splitText()
        then:
        result.val == 'foo\n'
        result.val == 'bar\n'
        result.val == Channel.STOP

        when:
        result = Channel.from('foo\nbar\nbaz').splitText(by:2)
        then:
        result.val == 'foo\nbar\n'
        result.val == 'baz\n'
        result.val == Channel.STOP
    }

    def 'should split text and invoke closure' () {

        when:
        def result = Channel.from('foo\nbar').splitText { it.trim().reverse() }
        then:
        result.val == 'oof'
        result.val == 'rab'
        result.val == Channel.STOP

        when:
        result = Channel.from('aa\nbb\ncc\ndd').splitText(by:2) { it.trim() }
        then:
        result.val == 'aa\nbb'
        result.val == 'cc\ndd'
        result.val == Channel.STOP
    }


}
