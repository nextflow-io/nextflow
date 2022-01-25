/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.extension

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class JoinOpTest extends Specification {

    def setup() {
        new Session()
    }

    def 'should join entries' () {
        given:
        def ch1 = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
        def ch2 = Channel.from(['Z', 6], ['Y', 5], ['X', 4])

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
    }



    def 'should join entries by index' () {
        given:
        def ch1 = Channel.from([1, 'X'], [2, 'Y'], [3, 'Z'], [7, 'P'])
        def ch2 = Channel.from([6, 'Z'], [5, 'Y'], [4, 'X'])

        when:
        def op = new JoinOp(ch1, ch2, [by:1])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
    }

    def 'should join entries with composite index' () {
        given:
        def ch1 = Channel.from([1, 'a','b', ['foo']], [2, 'p','q', ['bar']], [3, 'x','y', ['baz']], [7, 'P'])
        def ch2 = Channel.from([5, 'p','q', [333]], [4, 'a','b', [444]], [6, 'x','y', [555]])

        when:
        def op = new JoinOp(ch1, ch2, [by:[1,2]])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result.contains( ['a','b', 1, ['foo'], 4, [444]] )
        result.contains( ['p','q', 2, ['bar'], 5, [333]] )
        result.contains( ['x','y', 3, ['baz'], 6, [555]] )

    }


    def 'should join entries with remainder' () {
        given:
        def ch1 = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
        def ch2 = Channel.from(['Z', 6], ['Y', 5], ['X', 4], ['Q', ['foo','bar', [77,88,99]]])

        when:
        def op = new JoinOp(ch1, ch2, [remainder: true])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 5
        result.contains( ['X', 1, 4] )
        result.contains( ['Y', 2, 5] )
        result.contains( ['Z', 3, 6] )
        result.contains( ['P', 7, null] )
        result.contains( ['Q', null, ['foo','bar', [77,88,99]]])
    }

    def 'should join single item channels' () {

        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result.size() == 3
        result == [1,2,3]
    }

    def 'should join single item channels with remainder' () {

        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.from( 1,0,0,2,7,8,9,3 )

        when:
        def op = new JoinOp(ch1, ch2, [remainder: true])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 8
        result == [1, 2, 3, 0, 0, 7, 8, 9]
    }

    def 'should join empty channel and remainder' () {

        when:
        def left = Channel.from(1,2,3)
        def right = Channel.empty()
        def result = left.join(right, remainder: true)
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }

    def 'should join empty channel with pairs and remainder' () {

        when:
        def left = Channel.from(['X', 1], ['Y', 2], ['Z', 3])
        def right = Channel.empty()
        def result = left.join(right, remainder: true)
        then:
        result.val == ['X', 1, null]
        result.val == ['Y', 2, null]
        result.val == ['Z', 3, null]
        result.val == Channel.STOP
    }

    def 'should join a singleton value' () {

        when:
        given:
        def ch1 = Channel.from( 1,2,3 )
        def ch2 = Channel.value(1)

        when:
        def op = new JoinOp(ch1, ch2)
        def result = op.apply().toList().getVal()
        then:
        result == [1]
    }


    def 'should join pair with singleton and remainder' () {

        when:
        def left = Channel.from(['P', 0], ['X', 1], ['Y', 2], ['Z', 3])
        def right = Channel.from('X', 'Y', 'Z', 'Q')
        def result = left.join(right)
        then:
        result.val == ['X', 1]
        result.val == ['Y', 2]
        result.val == ['Z', 3]
        result.val == Channel.STOP

        when:
        left = Channel.from(['P', 0], ['X', 1], ['Y', 2], ['Z', 3])
        right = Channel.from('X', 'Y', 'Z', 'Q')
        result = left.join(right, remainder: true).toList().val.sort { it -> it[0] }
        then:
        result[2] == ['X', 1]
        result[3] == ['Y', 2]
        result[4] == ['Z', 3]
        result[0] == ['P', 0]
        result[1] == ['Q', null]

    }

    def 'should match gstrings' () {
        given:
        def A = "A"; def B = "B"; def C = "C"
        def left = Channel.of(['A', 'hola'], ['B', 'hello'], ['C', 'ciao'])
        def right = Channel.of(["$A", 'mundo'], ["$B", 'world'], ["$C", 'mondo'] )
        when:
        def result = left.join(right).toList().val.sort { it[0] }
        then:
        result[0] == ['A','hola','mundo']
        result[1] == ['B','hello','world']
        result[2] == ['C','ciao','mondo']

    }


    def 'should not fail on mismatches' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1], ['Y', 2])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 6], ['Y', 5])

        when:
        def op = new JoinOp(ch1, ch2, [failOnMismatch:true])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 2
        result.contains( ['X', 1, 6] )
        result.contains( ['Y', 2, 5] )
    }

    def 'should should fail on mismatches' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 6], ['Y', 5])
        and:
        def sess = Global.session as Session

        when:
        def op = new JoinOp(ch1, ch2, [failOnMismatch:true])
        def result = op.apply().toList().getVal()
        and:
        await(sess)
        then:
        sess.isAborted()
        sess.getError().message == 'Join mismatch for the following entries: key=Y values=[5]'
    }

    def 'should format error message '() {
        given:
        def op = new JoinOp(Mock(DataflowReadChannel), Mock(DataflowReadChannel), [:])

        when:
        op.checkForMismatch([:])
        then:
        noExceptionThrown()

        when:
        def buffer1 = [
                X: [:],
                Y: [(1): ['a','b']] ]
        op.checkForMismatch(buffer1)
        then:
        def e1 = thrown(AbortOperationException)
        e1.message == 'Join mismatch for the following entries: key=Y values=a,b'

        when:
        def buffer2 = [
                X: [(0): ['foo']],
                Y: [(1): ['a','b']] ]
        op.checkForMismatch(buffer2)
        then:
        def e2 = thrown(AbortOperationException)
        e2.message == 'Join mismatch for the following entries: \n- key=X values=foo \n- key=Y values=a,b'
    }

    def 'should not fail on duplicate matches' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1], ['X', 3])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 2], ['X', 4])

        when:
        def op = new JoinOp(ch1, ch2, [:])
        def result = op.apply().toList().getVal()
        then:
        result.size() == 2
        result.contains( ['X', 1, 2] )
        result.contains( ['X', 3, 4] )
    }

    def 'should fail on duplicate matches' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1], ['X', 3], ['X', 5])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 2], ['X', 4], ['X', 6])
        and:
        def sess = Global.session as Session

        when:
        def op = new JoinOp(ch1, ch2, [failOnDuplicate:true])
        def result = op.apply().toList().getVal()
        println "result=$result"
        and:
        await(sess)
        then:
        sess.isAborted()
        and:
        sess.error.message ==~ /Detected join operation duplicate emission on (left|right) channel -- offending element: key=X; value=(3|4|5|6)/
    }

    def 'should fail on duplicate with remainder' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1], ['X', 3])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 2])
        and:
        def sess = Global.session as Session

        when:
        def op = new JoinOp(ch1, ch2, [failOnDuplicate:true, remainder: true])
        def result = op.apply().toList().getVal()
        and:
        await(sess)
        then:
        sess.isAborted()
        sess.getError().message == 'Detected join operation duplicate emission on left channel -- offending element: key=X; value=3'
    }

    def 'should fail on duplicate without remainder' () {
        given:
        def ch1 = (DataflowReadChannel) Channel.of(['X', 1], ['X', 3])
        def ch2 = (DataflowReadChannel) Channel.of(['X', 2])
        and:
        def sess = Global.session as Session

        when:
        def op = new JoinOp(ch1, ch2, [failOnDuplicate:true])
        def result = op.apply().toList().getVal()
        then:
        await(sess)
        then:
        sess.isAborted()
        sess.getError().message == 'Detected join operation duplicate emission on left channel -- offending element: key=X; value=3'
    }


    protected void await(Session session) {
        def begin = System.currentTimeMillis()
        while( !session.isAborted() && System.currentTimeMillis()-begin<5_000 )
            sleep 100
    }
}
