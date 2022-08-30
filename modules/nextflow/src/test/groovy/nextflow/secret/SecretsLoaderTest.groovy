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

import spock.lang.Ignore
import spock.lang.Specification

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class SecretsLoaderTest extends Specification {

    def cleanupSpec() {
        def processEnvironmentClass = System.getenv().getClass()
        def field = processEnvironmentClass.getDeclaredField('m')
        field.accessible = true
        def map = (Map<String, String>) field.get(System.getenv())
        map.remove('NXF_ENABLE_SECRETS')
    }

    def "should check if NXF_ENABLE_SECRETS is #ENV"() {
        given:
        def processEnvironmentClass = System.getenv().getClass()
        def field = processEnvironmentClass.getDeclaredField('m')
        field.accessible = true
        def map = (Map<String, String>) field.get(System.getenv())
        if( ENV )
            map.put('NXF_ENABLE_SECRETS', "$ENV".toString())
        else
            map.remove('NXF_ENABLE_SECRETS')

        when:
        boolean enabled = SecretsLoader.isEnabled()
        then:
        enabled == RESULT

        where:
        ENV     | RESULT
        null    | true
        'true'  | true
        '1'     | false
        'valid' | false
        'false' | false
        'f'     | false
    }

}
