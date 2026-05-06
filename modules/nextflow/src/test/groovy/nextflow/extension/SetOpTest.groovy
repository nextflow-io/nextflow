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

import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class SetOpTest extends Dsl2Spec {

    def 'should set a channel in the global context' () {
        when:
        def result = runScript('''
            channel.of(1,2,3) | set { foo }
            foo | map { it * 2 }
        ''')
        then:
        result.val == 2
        result.val == 4
        result.val == 6

        when:
        result = runScript('''
            channel.value(5) | set { foo }
            foo | map { it * 2 }
        ''')
        then:
        result.val == 10
    }

    def 'should invoke set with dot notation' () {
        when:
        def result = runScript('''
            channel.of(1,2,3).set { foo }
            foo.map { it * 2 }
        ''')
        then:
        result.val == 2
        result.val == 4
        result.val == 6

        when:
        result = runScript('''
            channel.value('hello').set { foo }
            foo.map { it.toUpperCase() }
        ''')
        then:
        result.val == 'HELLO'
    }

}
