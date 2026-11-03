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
import nextflow.Global
import spock.lang.Specification

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowTapExtensionTest extends Specification {

    def 'should `tap` item to a new channel' () {

        when:
        def (result, first) = runDataflow {
            result = Channel.of( 4,7,9 ) .tap { first }.map { it+1 }
            [ result, first ]
        }
        then:
        first.val == 4
        first.val == 7
        first.val == 9
        first.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !Global.session.dag.isEmpty()

    }

    def 'should `tap` item to more than one channel' () {

        when:
        def (result, foo, bar) = runDataflow {
            result = Channel.of( 4,7,9 ) .tap { foo; bar }.map { it+1 }
            [ result, foo, bar ]
        }
        then:
        foo.val == 4
        foo.val == 7
        foo.val == 9
        foo.val == Channel.STOP
        bar.val == 4
        bar.val == 7
        bar.val == 9
        bar.val == Channel.STOP

        result.val == 5
        result.val == 8
        result.val == 10
        result.val == Channel.STOP

        !Global.session.dag.isEmpty()

    }

}
