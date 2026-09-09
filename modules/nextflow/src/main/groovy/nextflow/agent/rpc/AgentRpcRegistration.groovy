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

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.agent.AgentRunner

/** Connection material issued for one canonical agent task attempt. */
@Canonical
/*
 * @Canonical is a meta-annotation that implies @ToString over EVERY property, so the default render
 * puts the capability token in whatever consumed it -- a `log.debug "$registration"`, a Spock failure
 * message, an exception built by string interpolation. The broker goes out of its way never to log
 * the token, and one such line would undo that silently; an explicit @ToString wins over the one
 * @Canonical implies, so this is the whole fix. `includeNames` is not cosmetic either: without it the
 * render is four bare positional values with one of them missing, which reads as a malformed object
 * rather than a redacted one.
 */
@ToString(excludes='token', includeNames=true)
@CompileStatic
class AgentRpcRegistration {
    String invocationId
    String token
    String endpoint
    /**
     * SHA-256 of the broker's certificate, lowercase hex, which the task pins to authenticate the
     * driver; {@code null} when the broker serves cleartext (`agent.rpc.tls = false`). Unlike the
     * token this is a public commitment, not a secret, so putting it in the task script is harmless.
     */
    String fingerprint
    /**
     * Set only by a broker that deliberately serves cleartext (`agent.rpc.tls = false`). Cleartext is
     * carried as its own flag rather than inferred from a missing {@link #fingerprint} so that the
     * two cannot be confused: a runner that simply forgets the digest has a bug, not a licence to
     * dial unencrypted.
     */
    boolean insecure

    /**
     * The transport flags the proxy needs to dial the broker: pin the driver's certificate, or opt
     * out of transport security explicitly. An absent digest is never read as "unpinned" -- a
     * registration that carries neither fails here, on the driver, with a message naming the runner
     * contract, rather than sending the proxy to dial cleartext against a TLS listener and surfacing
     * as an unrelated connection failure inside the task.
     *
     * <p>When a registration somehow carries BOTH -- which the in-tree broker cannot produce, since it
     * derives {@link #insecure} from the same {@code agent.rpc.tls} switch that decides whether a
     * fingerprint exists at all -- the pin wins. Do not "fix" this into a throw for symmetry with the
     * proxy, which does reject the pair: there the two flags are its own argv and disagreeing about
     * them is a bug in the caller, whereas here the only tiebreak that cannot silently downgrade a
     * confused out-of-tree runner to cleartext is the secure one.
     */
    List<String> transportArgs() {
        if( fingerprint )
            return List.<String>of('--fingerprint', fingerprint)
        if( insecure )
            return List.<String>of('--insecure')
        throw new IllegalStateException("Agent RPC registration ${invocationId} carries no certificate fingerprint and does not opt out of transport security -- an AgentRunner serving TLS must return the broker's fingerprint".toString())
    }
}
