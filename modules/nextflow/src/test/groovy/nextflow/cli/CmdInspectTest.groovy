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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdInspectTest extends Specification {

    def 'should configure max rate' () {
        given:
        def cmd = new CmdInspect()

        when:
        def cfg1 = [:]
        cmd.configureMaxRate(cfg1)
        then:
        cfg1 == [wave:[httpClient:[maxRate:'5/30sec']]]

        when:
        def cfg2 = [wave:[enabled:true, httpClient: [something:true, maxRate: '1/s']]]
        cmd.configureMaxRate(cfg2)
        then:
        cfg2 == [wave:[enabled:true, httpClient: [something:true, maxRate: '5/30sec']]]
    }

}
