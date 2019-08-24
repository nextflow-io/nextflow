/*
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

import spock.lang.Timeout

import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class MixOp2Test extends Dsl2Spec {

    def 'should mix channels'() {
        when:
        def result = dsl_eval('''
            c1 = Channel.from( 1,2,3 )
            c2 = Channel.from( 'a','b' )
            c3 = Channel.value( 'z' )
            c1.mix(c2,c3)
            
        ''') .toList().val

        then:
        1 in result
        2 in result
        3 in result
        'a' in result
        'b' in result
        'z' in result
        !('c' in result)

    }

    def 'should mix with value channels'() {
        when:
        def result = dsl_eval('''
            Channel.value(1).mix( Channel.from([2,3])  )
            ''')
        then:
        result.toList().val.sort() == [1,2,3]
    }

    def 'should mix with two singleton'() {
        when:
        def result = dsl_eval('''
            Channel.value(1).mix( Channel.value(2)  )
            ''')
        then:
        result.toList().val.sort() == [1,2]
    }

}
