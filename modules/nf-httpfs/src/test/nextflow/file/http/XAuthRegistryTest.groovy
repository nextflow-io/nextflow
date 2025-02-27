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

package nextflow.file.http


import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class XAuthRegistryTest extends Specification {

    @Unroll
    def 'should register authenticator'() {
        given:
        def authenticator = new XAuthRegistry()
        and:
        def provider = Mock(XAuthProvider)
        authenticator.register(provider)
        and:
        def conn1 = Mock(URLConnection)

        when:
        def result = authenticator.authorize(conn1)
        then:
        1 * provider.authorize(conn1) >> ACCEPT
        and:
        result == EXPECTED

        cleanup:
        authenticator.unregister(provider)

        where:
        ACCEPT  | EXPECTED
        false   | false
        true    | true
    }
}
