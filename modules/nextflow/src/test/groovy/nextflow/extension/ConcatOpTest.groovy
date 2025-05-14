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
 */

package nextflow.extension

import spock.lang.Timeout

import nextflow.Channel
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ConcatOpTest extends Dsl2Spec {

    def 'should concat two channel'() {

        when:
        def result = dsl_eval('''
            c1 = Channel.of(1,2,3)
            c2 = Channel.of('a','b','c')
            c1.concat(c2)
        ''')
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == 3
        result.unwrap() == 'a'
        result.unwrap() == 'b'
        result.unwrap() == 'c'
        result.unwrap() == Channel.STOP
    }

    def 'should concat value with channel'() {
        when:
        def result = dsl_eval('''
            ch1 = Channel.value(1)
            ch2 = Channel.of(2,3)
            ch1.concat(ch2)        
        ''')
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == 3
        result.unwrap() == Channel.STOP
    }

    def 'should concat two value channels'() {
        when:
        def result = dsl_eval('''
            ch1 = Channel.value(1)
            ch2 = Channel.value(2)
            ch1.concat(ch2)        
        ''')
        then:
        result.unwrap() == 1
        result.unwrap() == 2
        result.unwrap() == Channel.STOP
    }

    def 'should concat with empty'() {
        when:
        def result = dsl_eval('''
            ch1 = Channel.value(1)
            ch2 = Channel.empty()
            ch1.concat(ch2)        
        ''')
        then:
        result.unwrap() == 1
        result.unwrap() == Channel.STOP
        
        when:
        result = dsl_eval('''
            ch1 = Channel.empty()
            ch2 = Channel.empty()
            ch1.concat(ch2)        
        ''')
        then:
        result.unwrap() == Channel.STOP
    }

}
