/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.extension

import nextflow.Channel
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class UntilManyOpTest extends Specification {

    def 'should emit channel items until the condition is verified' () {

        when:
        def source = Channel.of(1,2,3,4)
        def result = new UntilManyOp([source], { it==3 }).apply().get(0)
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == Channel.STOP

        when:
        source = Channel.of(1,2,3)
        result = new UntilManyOp([source], { it==5 }).apply().get(0)
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == 3
        result.unwrap() == Channel.STOP

    }

    def 'should emit channels util tuple condition is verified' () {
        given:
        def A = Channel.of(1,2,3)
        def B = Channel.of('alpha','beta')
        def C = Channel.of('foo', 'bar')

        when:
        def condition = { a, b, c -> a>2 }
        def (X,Y,Z) = new UntilManyOp([A,B,C], condition).apply()

        then:
        X.unwrap() == 1
        Y.unwrap() == 'alpha'
        Z.unwrap() == 'foo'
        and:
        X.unwrap() == 2
        Y.unwrap() == 'beta'
        Z.unwrap() == 'bar'
        and:
        X.unwrap() == Channel.STOP
        Y.unwrap() == Channel.STOP
        Z.unwrap() == Channel.STOP
    }

    def 'should emit channels until list condition is verified' () {
        given:
        def A = Channel.of(1,2,3)
        def B = Channel.of('alpha','beta')
        def C = Channel.of('foo', 'bar')

        when:
        def condition = { it -> it[0]>2 }
        def (X,Y,Z) = new UntilManyOp([A,B,C], condition).apply()

        then:
        X.unwrap() == 1
        Y.unwrap() == 'alpha'
        Z.unwrap() == 'foo'
        and:
        X.unwrap() == 2
        Y.unwrap() == 'beta'
        Z.unwrap() == 'bar'
        and:
        X.unwrap() == Channel.STOP
        Y.unwrap() == Channel.STOP
        Z.unwrap() == Channel.STOP
    }

}
