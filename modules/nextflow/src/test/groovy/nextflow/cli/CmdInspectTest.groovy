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
 *
 */

package nextflow.cli

import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdInspectTest extends Specification {

    def 'should ask for confirmation' () {
        given:
        def cmd = Spy(new CmdInspect(awaitMode: AWAIT))
        Map wave

        when:
        wave = CONFIG
        cmd.checkWaveConfig(wave)
        then:
        INVOCTIONS * cmd.promptConfirmation() >> REPLY
        and:
        wave == EXPECTED

        where:
        CONFIG                              | INVOCTIONS    | REPLY     | AWAIT     | EXPECTED
        [:]                                 | 0             | null      | false     | [:]
        [enabled: true]                     | 0             | null      | false     | [enabled: true]
        [enabled: true, freeze: true]       | 1             | 'Y'       | false     | [enabled: true, freeze: true, awaitMode: false]
        [enabled: true, freeze: true]       | 1             | 'Y'       | true      | [enabled: true, freeze: true, awaitMode: true]

    }

    def 'should abort the operation' () {
        given:
        def cmd = Spy(CmdInspect)

        when:
        cmd.checkWaveConfig([enabled:true,freeze: true])
        then:
        1 * cmd.promptConfirmation() >> 'n'
        and:
        thrown(AbortOperationException)
    }
}
