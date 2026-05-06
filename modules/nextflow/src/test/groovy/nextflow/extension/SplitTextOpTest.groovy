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
class SplitTextOpTest extends Specification {

    def 'should split text' () {

        when:
        def result = runDataflow {
            Channel.of('foo\nbar').splitText()
        }
        then:
        result.val == 'foo\n'
        result.val == 'bar\n'
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of('foo\nbar\nbaz').splitText(by:2)
        }
        then:
        result.val == 'foo\nbar\n'
        result.val == 'baz\n'
        result.val == Channel.STOP
    }

    def 'should split text and invoke closure' () {

        when:
        def result = runDataflow {
            Channel.of('foo\nbar').splitText { it.trim().reverse() }
        }
        then:
        result.val == 'oof'
        result.val == 'rab'
        result.val == Channel.STOP

        when:
        result = runDataflow {
            Channel.of('aa\nbb\ncc\ndd').splitText(by:2) { it.trim() }
        }
        then:
        result.val == 'aa\nbb'
        result.val == 'cc\ndd'
        result.val == Channel.STOP
    }


}
