/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowMergeExtensionTest extends Specification {

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }

    def 'should merge with open array with custom closure'() {
        when:
        def alpha = Channel.from(1, 3, 5);
        def beta =  Channel.from(2, 4, 6);
        def delta = Channel.from(7, 8, 1);
        def result = alpha.merge( beta, delta ) { a,b,c -> [c,b,a] }
        then:
        result instanceof DataflowQueue
        result.val == [7,2,1]
        result.val == [8,4,3]
        result.val == [1,6,5]
        result.val == Channel.STOP
    }

    def 'should merge with open array' () {
        when:
        def alpha = Channel.from(1, 3, 5);
        def beta =  Channel.from(2, 4, 6);
        def delta = Channel.from(7, 8, 1);
        def result = alpha.merge( beta, delta )
        then:
        result instanceof DataflowQueue
        result.val == [1,2,7]
        result.val == [3,4,8]
        result.val == [5,6,1]
        result.val == Channel.STOP
    }

    def 'should merge with with default'() {

        when:
        def left =  Channel.from(1, 3, 5);
        def right = Channel.from(2, 4, 6);
        def result = left.merge(right)
        then:
        result instanceof DataflowQueue
        result.val == [1,2]
        result.val == [3,4]
        result.val == [5,6]
        result.val == Channel.STOP

        when:
        left  = Channel.from(1, 2, 3);
        right = Channel.from(['a','b'], ['p','q'], ['x','z']);
        result = left.merge(right)
        then:
        result instanceof DataflowQueue
        result.val == [1, 'a','b']
        result.val == [2, 'p','q']
        result.val == [3, 'x','z']
        result.val == Channel.STOP

        when:
        left  = Channel.from('A','B','C');
        right = Channel.from(['a',[1,2,3]], ['b',[3,4,5]], ['c',[6,7,8]]);
        result = left.merge(right)
        then:
        result instanceof DataflowQueue
        result.val == ['A', 'a', [1,2,3]]
        result.val == ['B', 'b', [3,4,5]]
        result.val == ['C', 'c', [6,7,8]]
        result.val == Channel.STOP

    }

    def 'should merge with list'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta  = Channel.from(2, 4, 6);
        def delta = Channel.from(7, 8, 1);
        def result = alpha.merge( [beta, delta] ) { a,b,c -> [c,b,a] }
        then:
        result instanceof DataflowQueue
        result.val == [7,2,1]
        result.val == [8,4,3]
        result.val == [1,6,5]
        result.val == Channel.STOP

    }

    def 'should merge with queue'() {

        when:
        def alpha = Channel.from(1, 3, 5);
        def beta = Channel.from(2, 4, 6);

        def result = alpha.merge(beta) { a,b -> [a, b+1] }

        then:
        result instanceof DataflowQueue
        result.val == [1,3]
        result.val == [3,5]
        result.val == [5,7]
        result.val == Channel.STOP
    }

    def 'should merge with variables with custom closure'() {

        when:
        def alpha = Channel.value('Hello');
        def beta = Channel.value('World')
        def result = alpha.merge(beta) { a,b -> [b, a] }
        then:
        result instanceof DataflowVariable
        result.val == ['World', 'Hello']
        result.val == ['World', 'Hello']
    }

    def 'should merge variables' () {
        when:
        def alpha = Channel.value('Hello');
        def beta = Channel.value('World')
        def result = alpha.merge(beta)
        then:
        result instanceof DataflowVariable
        result.val == ['Hello','World']
        result.val == ['Hello','World']
    }

}
