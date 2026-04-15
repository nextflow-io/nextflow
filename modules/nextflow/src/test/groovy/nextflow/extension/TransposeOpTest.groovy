/*
 * Copyright 2013-2026, Seqera Labs
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

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import spock.lang.Specification
import spock.lang.Timeout

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class TransposeOpTest extends Specification {

    def 'should transpose tuple' () {

        given:
        def values = [ ['a',[1,2,3],'p','q'], ['b',[4,5,6],'x','y'] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose()
        }
        then:
        result.val == ['a',1,'p','q']
        result.val == ['a',2,'p','q']
        result.val == ['a',3,'p','q']

        result.val == ['b',4,'x','y']
        result.val == ['b',5,'x','y']
        result.val == ['b',6,'x','y']

        result.val == Channel.STOP
    }

    def 'should transpose multiple tuples' () {

        given:
        def values = [ ['a',[1,2,3],['p','q']], ['b',[4,5,6],['x','y']] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose()
        }
        then:
        result.val == ['a',1,'p']
        result.val == ['a',2,'q']

        result.val == ['b',4,'x']
        result.val == ['b',5,'y']

        result.val == Channel.STOP
    }

    def 'should transpose multiple tuples with remainder' () {

        given:
        def values = [ ['a',[1,2,3],['p','q']], ['b',[4,5],['x','y','z']] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose(remainder: true)
        }
        then:
        result.val == ['a',1,'p']
        result.val == ['a',2,'q']
        result.val == ['a',3,null]

        result.val == ['b',4,'x']
        result.val == ['b',5,'y']
        result.val == ['b',null,'z']

        result.val == Channel.STOP
    }

    def 'should transpose tuple by 1' () {

        given:
        def values = [ ['a',[1,2,3],['p','q']], ['b',[4,5,6],['x','y']] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose(by: 1)
        }
        then:
        result.val == ['a',1,['p','q']]
        result.val == ['a',2,['p','q']]
        result.val == ['a',3,['p','q']]

        result.val == ['b',4,['x','y']]
        result.val == ['b',5,['x','y']]
        result.val == ['b',6,['x','y']]

        result.val == Channel.STOP
    }

    def 'should transpose tuple with index' () {

        given:
        def values = [ ['a',[1,2,3],['p','q']], ['b',[4,5,6],['x','y']] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose(by: 2)
        }
        then:
        result.val == ['a',[1,2,3],'p']
        result.val == ['a',[1,2,3],'q']
        result.val == ['b',[4,5,6],'x']
        result.val == ['b',[4,5,6],'y']
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.fromList(values).transpose(by: 1)
        }
        then:
        result.val == ['a',1,['p','q']]
        result.val == ['a',2,['p','q']]
        result.val == ['a',3,['p','q']]
        result.val == ['b',4,['x','y']]
        result.val == ['b',5,['x','y']]
        result.val == ['b',6,['x','y']]
        result.val == Channel.STOP

    }

    def 'should transpose tuple 3' () {

        given:
        def values = [ ['a','b'], ['c','d'], ['e','f'] ]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose()
        }
        then:
        result.val == ['a','b']
        result.val == ['c', 'd']
        result.val == ['e', 'f']

        result.val == Channel.STOP
    }

    def 'should transpose values' () {

        given:
        def values = ['a','b','c','d']

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose()
        }
        then:
        result.val == 'a'
        result.val == 'b'
        result.val == 'c'
        result.val == 'd'
        result.val == Channel.STOP
    }

    def 'should transpose a tuple' () {
        given:
        def target = Mock(DataflowWriteChannel)
        def op = new TransposeOp(target:target)

        when:
        op.transpose0(['a','b','c'])
        then:
        1 * target.bind(['a','b','c'])

        when:
        op.transpose0(['a','b',['p','q']])
        then:
        1 * target.bind(['a','b','p'])
        then:
        1 * target.bind(['a','b','q'])

        when:
        op.transpose0(['a', 'b', ['p','q'], 'c', ['x','y','z']])
        then:
        1 * target.bind(['a','b','p','c','x'])
        then:
        1 * target.bind(['a','b','q','c','y'])

        when:
        op.transpose0(['a','b',[], ['p','q']])
        then:
        0 * target.bind(_)

    }

    def 'should raise an exception' () {
        given:
        def target = Mock(DataflowWriteChannel)
        def op = new TransposeOp(target:target, cols:[0])

        when:
        op.transpose0(['a',['b','c']])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should transpose a tuple with reminder' () {
        given:
        def target = Mock(DataflowWriteChannel)
        def op = new TransposeOp(target:target, remainder: true)

        when:
        op.transpose0(['a', 'b', ['p','q'], 'c', ['x','y','z']])
        then:
        1 * target.bind(['a','b', 'p', 'c','x'])
        then:
        1 * target.bind(['a','b', 'q', 'c','y'])
        then:
        1 * target.bind(['a','b', null,'c','z'])

        when:
        op.transpose0(['a', 'b', ['p'], 'c', ['x','y','z']])
        then:
        1 * target.bind(['a','b', 'p',  'c', 'x'])
        then:
        1 * target.bind(['a','b', null, 'c', 'y'])
        then:
        1 * target.bind(['a','b', null, 'c', 'z'])

        when:
        op.transpose0(['a', 'b', [], 'c', ['x','y','z']])
        then:
        1 * target.bind(['a','b', null, 'c', 'x'])
        then:
        1 * target.bind(['a','b', null, 'c', 'y'])
        then:
        1 * target.bind(['a','b', null, 'c', 'z'])
    }

    def 'should fetch indexes from a tuple' () {
        given:
        def ch = Mock(DataflowQueue)
        def op = new TransposeOp(ch)

        expect:
        op.findIndexes(['a','b','c']) == []
        op.findIndexes(['a','b',['c','d']]) == [2]
        op.findIndexes(['a',['b'],'x',['c','d']]) == [1,3]
    }

    def 'should apply transpose operator' () {

        given:
        def values = [['a',[1,2],['x']], ['b',[3,4], ['p','q']]]

        when:
        def result = runDataflow {
            Channel.fromList(values).transpose()
        }
        then:
        result.val == ['a',1,'x']
        result.val == ['b',3,'p']
        result.val == ['b',4,'q']
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.fromList(values).transpose(by:1)
        }
        then:
        result.val == ['a',1,['x']]
        result.val == ['a',2,['x']]
        result.val == ['b',3,['p','q']]
        result.val == ['b',4,['p','q']]
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.fromList(values).transpose(remainder:true)
        }
        then:
        result.val == ['a',1,'x']
        result.val == ['a',2,null]
        result.val == ['b',3,'p']
        result.val == ['b',4,'q']
        result.val == Channel.STOP
    }
}
