/*
 * Copyright 2020-2022, Seqera Labs
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

package io.seqera.wave.plugin

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveFactoryTest extends Specification {

    def 'should not change config/1' () {
        given:
        def session = Mock(Session) { getConfig() >> CONFIG }
        def factory = new WaveFactory()

        when:
        factory.create(session)
        then:
        CONFIG == EXPECTED

        where:
        CONFIG                      | EXPECTED
        [:]                         | [:]
        [wave:[enabled:true]]       | [wave:[enabled:true]]
        [wave:[enabled:true], fusion:[enabled:true]]     | [disableRemoteBinDir: true, wave:[enabled:true,bundleProjectResources:true], fusion:[enabled:true]]
    }

}
