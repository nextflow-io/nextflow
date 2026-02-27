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

import nextflow.Channel
import spock.lang.Specification

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class UntilManyOpTest extends Specification {

    def 'should emit channel items until the condition is verified' () {

        when:
        def result = runDataflow {
            def source = Channel.of(1,2,3,4)
            new UntilManyOp([read(source)], { it==3 }).apply().get(0)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == Channel.STOP

        when:
        result = runDataflow {
            def source = Channel.of(1,2,3)
            new UntilManyOp([read(source)], { it==5 }).apply().get(0)
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP

    }

    def 'should emit channels util tuple condition is verified' () {

        when:
        def (X,Y,Z) = runDataflow {
            def A = Channel.of(1,2,3)
            def B = Channel.of('alpha','beta')
            def C = Channel.of('foo', 'bar')
            def condition = { a, b, c -> a>2 }
            new UntilManyOp([read(A), read(B), read(C)], condition).apply()
        }

        then:
        X.val == 1
        Y.val == 'alpha'
        Z.val == 'foo'
        and:
        X.val == 2
        Y.val == 'beta'
        Z.val == 'bar'
        and:
        X.val == Channel.STOP
        Y.val == Channel.STOP
        Z.val == Channel.STOP
    }

    def 'should emit channels until list condition is verified' () {

        when:
        def (X,Y,Z) = runDataflow {
            def A = Channel.of(1,2,3)
            def B = Channel.of('alpha','beta')
            def C = Channel.of('foo', 'bar')
            def condition = { abc -> abc[0]>2 }
            new UntilManyOp([read(A), read(B), read(C)], condition).apply()
        }

        then:
        X.val == 1
        Y.val == 'alpha'
        Z.val == 'foo'
        and:
        X.val == 2
        Y.val == 'beta'
        Z.val == 'bar'
        and:
        X.val == Channel.STOP
        Y.val == Channel.STOP
        Z.val == Channel.STOP
    }

    def read(channel) {
        channel.createReadChannel()
    }

}
