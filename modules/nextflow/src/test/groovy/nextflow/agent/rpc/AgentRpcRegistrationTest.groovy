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
package nextflow.agent.rpc

import spock.lang.Specification
import nextflow.agent.AgentRunner

/**
 * Pins the two properties of a registration that are security decisions rather than plumbing:
 * a missing certificate fingerprint must never be read as "dial cleartext", and the capability
 * token must never reach a log line, an exception or a test report through {@code toString()}.
 */
class AgentRpcRegistrationTest extends Specification {

    def 'should pin the driver certificate when a fingerprint is issued'() {
        given:
        def registration = new AgentRpcRegistration('inv-1', 'tok-1', '127.0.0.1:9999', 'abc123')

        expect:
        registration.transportArgs() == ['--fingerprint', 'abc123']
    }

    def 'should opt out of transport security only when the broker says so explicitly'() {
        given: 'the shape a broker built with `agent.rpc.tls = false` returns'
        def registration = new AgentRpcRegistration('inv-1', 'tok-1', '127.0.0.1:9999', null, true)

        expect:
        registration.transportArgs() == ['--insecure']
    }

    def 'should fail closed when a runner returns neither a pin nor an explicit opt-out'() {
        given: 'an AgentRunner that forgot the fingerprint -- register() has a default impl and this type is @Canonical, so the short constructor still compiles'
        def registration = new AgentRpcRegistration('inv-1', 'tok-1', '127.0.0.1:9999')

        when:
        registration.transportArgs()

        then: 'the driver refuses to build the command, rather than telling the proxy to dial cleartext at a TLS listener'
        // Downgrading silently would surface inside the task as an unrelated "connect to driver"
        // failure -- the same class of misleading error transport security exists to remove.
        def error = thrown(IllegalStateException)
        error.message.contains('inv-1')
        error.message.contains('fingerprint')
    }

    def 'should keep the capability token out of the rendered form'() {
        given:
        def registration = new AgentRpcRegistration('inv-1', 'super-secret-token', '127.0.0.1:9999', 'abc123')

        when:
        def rendered = registration.toString()

        then: 'the token is absent -- @Canonical implies a @ToString over EVERY property, so without the explicit one this renders the secret'
        !rendered.contains('super-secret-token')

        and: 'the remaining values are NAMED, so a redacted render is not mistaken for a malformed one'
        // Without `includeNames` the output is four bare positional values with one silently
        // missing -- which reads as a bug in the object, not as a deliberate redaction. This clause
        // is what pins that half: `@ToString(excludes='token')` alone still passes every other
        // assertion here.
        rendered.contains('invocationId:inv-1')
        rendered.contains('endpoint:127.0.0.1:9999')
        rendered.contains('fingerprint:abc123')
        rendered.contains('insecure:false')

        and: 'while everything else @Canonical generates is untouched'
        registration == new AgentRpcRegistration('inv-1', 'super-secret-token', '127.0.0.1:9999', 'abc123')
        registration.hashCode() == new AgentRpcRegistration('inv-1', 'super-secret-token', '127.0.0.1:9999', 'abc123').hashCode()
    }
}
