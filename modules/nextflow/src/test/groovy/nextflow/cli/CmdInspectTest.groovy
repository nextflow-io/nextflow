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

package nextflow.cli

import nextflow.container.inspect.ContainerInspectMode
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
        cmd.setupInspectMode()
        cmd.checkWaveDryRun(wave)
        def expectedDryRun = ContainerInspectMode.waveDryRun()
        cmd.disableInspectMode()
        then:
        expectedDryRun == EXPECTED_DRYRUN

        where:
        WAVE                            | CONCRETIZE    | EXPECTED_DRYRUN
        [:]                             | false         | true
        [:]                             | true          | true
        and:
        [enabled:true]                  | false         | true
        [enabled:true]                  | true          | true
        and:
        [enabled:true, freeze: true]    | false         | true
        [enabled:true, freeze: true]    | true          | false

    }

}
