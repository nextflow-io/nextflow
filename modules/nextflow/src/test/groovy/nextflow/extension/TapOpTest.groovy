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

import spock.lang.Timeout
import test.Dsl2Spec

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(10)
class TapOpTest extends Dsl2Spec {

    def 'should forward a queue channel into multiple targets' () {
        when:
        def result = dsl_eval('''
            result = Channel.of( 4,7,9 ).tap { foo; bar }.map { v -> v + 1 }
            [ foo, bar, result ]
            ''')

        then:
        result[0].getVal() == 4
        result[0].getVal() == 7
        result[0].getVal() == 9
        result[0].getVal() == Channel.STOP
        result[1].getVal() == 4
        result[1].getVal() == 7
        result[1].getVal() == 9
        result[1].getVal() == Channel.STOP
        result[2].getVal() == 5
        result[2].getVal() == 8
        result[2].getVal() == 10
        result[2].getVal() == Channel.STOP
    }

    def 'should forward a value channel into multiple targets' () {
        when:
        def result = dsl_eval('''
            Channel.value(4).tap { foo ; bar }
            [ foo, bar ]
            ''')
        then:
        result[0].getVal() == 4
        result[1].getVal() == 4
    }

}
