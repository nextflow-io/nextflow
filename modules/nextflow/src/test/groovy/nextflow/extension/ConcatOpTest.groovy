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

import nextflow.Channel
import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ConcatOpTest extends Dsl2Spec {

    def 'should concat two channel'() {

        when:
        def result = runScript('''
            c1 = channel.of(1,2,3)
            c2 = channel.of('a','b','c')
            c1.concat(c2)
        ''')
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == 'a'
        result.val == 'b'
        result.val == 'c'
        result.val == Channel.STOP
    }

    def 'should concat value with channel'() {
        when:
        def result = runScript('''
            ch1 = channel.value(1)
            ch2 = channel.of(2,3)
            ch1.concat(ch2)
        ''')
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
    }

    def 'should concat two value channels'() {
        when:
        def result = runScript('''
            ch1 = channel.value(1)
            ch2 = channel.value(2)
            ch1.concat(ch2)
        ''')
        then:
        result.val == 1
        result.val == 2
        result.val == Channel.STOP
    }

    def 'should concat with empty'() {
        when:
        def result = runScript('''
            ch1 = channel.value(1)
            ch2 = channel.empty()
            ch1.concat(ch2)
        ''')
        then:
        result.val == 1
        result.val == Channel.STOP

        when:
        result = runScript('''
            ch1 = channel.empty()
            ch2 = channel.empty()
            ch1.concat(ch2)
        ''')
        then:
        result.val == Channel.STOP
    }

}
