/*
 * Copyright 2023, Seqera Labs
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

package nextflow.script.params

import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ArityParamTest extends Specification {

    static class DefaultArityParam implements ArityParam {
        DefaultArityParam() {}
    }

    @Unroll
    def testArity () {

        when:
        def param = new DefaultArityParam()
        param.setArity(VALUE)
        then:
        param.arity.min == MIN
        param.arity.max == MAX
        param.isSingle() == SINGLE

        where:
        VALUE  | SINGLE | MIN | MAX
        '1'    | true   | 1   | 1
        '1..*' | false  | 1   | Integer.MAX_VALUE
        '0..*' | false  | 0   | Integer.MAX_VALUE
    }

    @Unroll
    def testArityRange () {

        when:
        def range = new ArityParam.Range(MIN, MAX)
        then:
        range.contains(2) == TWO
        range.toString() == STRING

        where:
        MIN | MAX               | TWO   | STRING
        1   | 1                 | false | '1'
        1   | Integer.MAX_VALUE | true  | '1..*'
    }

}
