/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class SecretsLoaderTest extends Specification {

    def "should be enabled by default" () {
        given:
        SysEnv.push([:])
        expect:
        SecretsLoader.isEnabled()
        cleanup:
        SysEnv.pop()
    }

    @Unroll
    def "should check if NXF_ENABLE_SECRETS is #ENV"() {
        given:
        SysEnv.push(NXF_ENABLE_SECRETS: "$ENV")

        when:
        boolean enabled = SecretsLoader.isEnabled()
        then:
        enabled == RESULT

        cleanup:
        SysEnv.pop()

        where:
        ENV     | RESULT
        'true'  | true
        '1'     | false
        'valid' | false
        'false' | false
        'f'     | false
    }

}
