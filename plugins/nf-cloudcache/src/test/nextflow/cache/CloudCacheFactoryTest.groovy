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

package nextflow.cache

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CloudCacheFactoryTest extends Specification {

    @Unroll
    def 'should be enabled only when the cloud cache is configured: #CONFIG'() {
        expect:
        new CloudCacheFactory().isEnabled(CONFIG) == EXPECTED

        where:
        CONFIG                                  | EXPECTED
        [cloudcache: [enabled: true]]           | true
        [cloudcache: [enabled: false]]          | false
        [cloudcache: [path: 's3://foo']]        | false
        [:]                                     | false
        null                                    | false
    }

}
