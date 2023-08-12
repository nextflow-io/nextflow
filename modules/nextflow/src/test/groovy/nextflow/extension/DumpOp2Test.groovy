/*
 * Copyright 2013-2023, Seqera Labs
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
import nextflow.Session
import test.Dsl2Spec
import test.OutputCapture

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class DumpOp2Test extends Dsl2Spec {

    @org.junit.Rule
    OutputCapture capture = new OutputCapture()

    def 'should support pipe' () {

        given:
        new Session(dumpChannels: ['*'])

        when:
        def result = dsl_eval(/
            Channel.of(1, 2, 3) | dump
        /)
        def stdout = capture.toString()
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
        stdout.contains('[DUMP] 1')
        stdout.contains('[DUMP] 2')
        stdout.contains('[DUMP] 3')

    }

}
