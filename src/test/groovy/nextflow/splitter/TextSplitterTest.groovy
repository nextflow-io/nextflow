package nextflow.splitter

import groovyx.gpars.dataflow.operator.PoisonPill
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TextSplitterTest extends Specification {

    def testTextByLine() {

        expect:
        new TextSplitter().target("Hello\nworld\n!").list() == ['Hello\n','world\n','!\n']
        new TextSplitter().target("Hello\nworld\n!").options(each:{ it.trim().reverse() }) .list()  == ['olleH','dlrow','!']

    }


    def testSplitLinesByCount () {

        expect:
        new TextSplitter().target("Hello\nHola\nHalo").list() == ['Hello\n', 'Hola\n', 'Halo\n']
        new TextSplitter().options(by:3).target("11\n22\n33\n44\n55").list() == [ '11\n22\n33\n', '44\n55\n' ]
        new TextSplitter().options(by:2).target("11\n22\n33\n44\n55").list() == [ '11\n22\n', '33\n44\n', '55\n' ]
        new TextSplitter().options(by:2).target('Hello\nworld\n!').list() == ['Hello\nworld\n','!\n']
    }

    def testSplitFileByLine () {

        setup:
        def file = File.createTempFile('chunk','test')
        file.deleteOnExit()
        file.text = '''\
        line1
        line2
        line3
        line4
        line5
        '''.stripIndent()

        when:
        def lines = new TextSplitter().target(file).list()

        then:
        lines[0] == 'line1\n'
        lines[1] == 'line2\n'
        lines[2]== 'line3\n'
        lines[3] == 'line4\n'
        lines[4] == 'line5\n'

        when:
        def channel = new TextSplitter().target(file).options(by:2).channel()
        then:
        channel.val == 'line1\nline2\n'
        channel.val == 'line3\nline4\n'
        channel.val == 'line5\n'

    }

    def testSplitChannel() {


        when:   '*into* is a Dataflow channel, get the chopped items, closing by a poison-pill'
        def channel = new TextSplitter().target("Hello\nworld\n!").channel()
        then:
        channel.val == 'Hello\n'
        channel.val == 'world\n'
        channel.val == '!\n'
        channel.val == PoisonPill.instance

    }
}
