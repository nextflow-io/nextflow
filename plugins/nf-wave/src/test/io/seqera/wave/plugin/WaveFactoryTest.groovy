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

package io.seqera.wave.plugin

import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveFactoryTest extends Specification {

    @Unroll
    def 'should not change config' () {
        given:
        def session = Mock(Session) { getConfig() >> CONFIG }
        def factory = new WaveFactory()

        when:
        factory.create(session)
        then:
        CONFIG == EXPECTED
        and:
        DISABLED * session.setDisableRemoteBinDir(true) >> null

        where:
        CONFIG                      | EXPECTED                  | DISABLED
        [:]                         | [:]                       | 0
        [wave:[enabled:true]]       | [wave:[enabled:true]]     | 0
        [wave:[enabled:true], fusion:[enabled:true]]     | [wave:[enabled:true,bundleProjectResources:true], fusion:[enabled:true]] | 1
    }

    @Unroll
    def 'should check s5cmd is enabled' () {
        given:
        def factory = new WaveFactory()

        expect:
        factory.isAwsBatchFargateMode(CONFIG) == EXPECTED
        
        where:
        CONFIG                                  | EXPECTED
        [:]                                     | false
        [aws:[batch:[platformType:'foo']]]       | false
        [aws:[batch:[platformType:'fargate']]]   | true
        [aws:[batch:[platformType:'Fargate']]]   | true

    }
    def 'should fail when wave is disabled' () {
        given:
        def CONFIG = [wave:[:], fusion:[enabled:true]]
        def session = Mock(Session) { getConfig() >> CONFIG }
        def factory = new WaveFactory()

        when:
        factory.create(session)
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Fusion feature requires enabling Wave service'
    }

    def 'should not fail when wave is disabled' () {
        given:
        SysEnv.push(NXF_DISABLE_WAVE_SERVICE: 'true')
        def CONFIG = [wave:[:], fusion:[enabled:true]]
        def session = Mock(Session) { getConfig() >> CONFIG }
        def factory = new WaveFactory()

        when:
        factory.create(session)
        then:
        noExceptionThrown()
        and:
        0 * session.setDisableRemoteBinDir(true) >> null
        and:
        CONFIG == [wave:[:], fusion:[enabled:true]]
    }


}
