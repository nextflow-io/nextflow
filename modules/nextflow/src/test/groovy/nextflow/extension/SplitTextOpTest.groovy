package nextflow.extension

import nextflow.Session
import spock.lang.Specification

import nextflow.Channel

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitTextOpTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should split text' () {
        when:
        def result = Channel.of('foo\nbar').splitText()
        then:
        result.unwrap() == 'foo\n'
        result.unwrap() == 'bar\n'
        result.unwrap() == Channel.STOP

        when:
        result = Channel.of('foo\nbar\nbaz').splitText(by:2)
        then:
        result.unwrap() == 'foo\nbar\n'
        result.unwrap() == 'baz\n'
        result.unwrap() == Channel.STOP
    }

    def 'should split text and invoke closure' () {
        when:
        def result = Channel.of('foo\nbar').splitText { it.trim().reverse() }
        then:
        result.unwrap() == 'oof'
        result.unwrap() == 'rab'
        result.unwrap() == Channel.STOP

        when:
        result = Channel.of('aa\nbb\ncc\ndd').splitText(by:2) { it.trim() }
        then:
        result.unwrap() == 'aa\nbb'
        result.unwrap() == 'cc\ndd'
        result.unwrap() == Channel.STOP
    }

}
