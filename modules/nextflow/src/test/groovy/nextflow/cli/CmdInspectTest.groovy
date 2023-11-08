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


import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdInspectTest extends Specification {

    @Unroll
    def 'should ask for confirmation' () {
        given:
        def cmd = Spy(new CmdInspect(concretize: CONCRETIZE))
        Map wave

        when:
        wave = WAVE
        cmd.checkWaveConfig(wave)
        then:
        wave == EXPECTED

        where:
        WAVE                            | CONCRETIZE    | EXPECTED
        [:]                             | false         | [:]
        [:]                             | true          | [:]
        and:
        [enabled:true]                  | false         | [enabled:true]
        [enabled:true]                  | true          | [enabled:true]
        and:
        [enabled:true, freeze: true]    | false         | [enabled:true, freeze:true, dryRun: true]
        [enabled:true, freeze: true]    | true          | [enabled:true, freeze:true, dryRun: false]

    }

}
