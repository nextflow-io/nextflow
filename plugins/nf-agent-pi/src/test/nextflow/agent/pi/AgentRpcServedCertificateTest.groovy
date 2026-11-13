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
package nextflow.agent.pi

import java.security.MessageDigest
import java.security.cert.Certificate
import java.security.cert.X509Certificate
import javax.net.ssl.SSLContext
import javax.net.ssl.SSLSocket
import javax.net.ssl.TrustManager
import javax.net.ssl.X509TrustManager

import nextflow.Global
import nextflow.Session
import nextflow.agent.AgentRunnerRequest
import spock.lang.Specification

/**
 * Closes the one loop that neither transport-security suite closes: that the certificate DER the
 * server actually puts ON THE WIRE is byte-identical to the DER the advertised fingerprint was
 * computed over.
 *
 * <p>Both existing suites compare against a same-language counterpart, which is why the gap
 * survives them. {@code agent-rpc}'s {@code main_test.go} builds its OWN Go certificate and stands
 * up a Go gRPC server, so it is a hand-maintained mirror of what {@link AgentRpcTlsCredentials}
 * emits, and mirrors drift. {@code AgentRpcBrokerTest} does dial the real broker, but with a
 * grpc-java client that installs the broker's PEM as a TRUST ANCHOR -- which validates by path
 * building, so a certificate that was re-encoded on the way out yet stayed semantically identical
 * would still be accepted there while hashing to a different digest.
 *
 * <p>That distinction is load-bearing because the driver never hands the task a certificate, only a
 * digest of one. grpc-netty-shaded bundles netty-tcnative, so on the platforms it ships natives for
 * -- CI included -- the PEM is re-parsed and re-serialized by BoringSSL before it reaches the
 * socket, rather than the PEM body being echoed through. Reproducing identical bytes for DER input
 * is what any conforming encoder does, so the risk is low; but if it ever were not so, BOTH existing
 * suites stay green while every real agent task dies with
 * {@code driver TLS certificate fingerprint mismatch} -- a message that reads like an attack and
 * would be diagnosed like one.
 *
 * <p>The probe is a bare JSSE {@link SSLSocket} rather than a shell-out to {@code openssl s_client}:
 * it always runs, so there is no skip to go stale; it parses no command output; and it introduces no
 * second TLS implementation as a variable. The assertion is over {@code chain[0].encoded} from the
 * PEER-RECEIVED certificate, i.e. bytes that came back from the server. Hashing
 * {@code broker.certificatePem} instead would restate what the broker already believes and prove
 * nothing new, which is why that accessor is deliberately untouched here.
 *
 * <p>What this does NOT close: the socket negotiates ALPN and then hangs up without ever sending an
 * HTTP/2 preface, so no h2 framing and no gRPC call is exercised. The cross-LANGUAGE half of the
 * handshake stays covered only by agent-rpc's Go tests against a Go server, and the h2-over-TLS half
 * by {@code AgentRpcBrokerTest}'s grpc-java client. This file is the byte-identity half alone;
 * retiring either of the others on the strength of it would lose real coverage.
 */
class AgentRpcServedCertificateTest extends Specification {

    /**
     * How long the handshake may take before the spec fails instead of hanging the suite. Generous
     * because it covers a loopback connect and one EC handshake, nothing more.
     */
    private static final int HANDSHAKE_TIMEOUT_MILLIS = 20_000

    AgentRpcBroker broker

    def setup() {
        Global.session = new Session([:])
        broker = AgentRpcBroker.get()
    }

    def cleanup() {
        broker?.close()
        Global.session = null
    }

    def 'should serve on the wire the exact certificate bytes the advertised fingerprint pins'() {
        given: 'a capability registered against the default broker, i.e. transport security on'
        def registration = broker.register(
            new AgentRunnerRequest(model: 'openai/test', prompt: 'confidential', requestTimeoutSeconds: 30), false)

        and: 'the broker really did mint a pin'
        // Diagnostic, not the property under test, and an EXPLICIT assert because a bare condition in
        // a setup block is only an expression -- Spock asserts implicitly in expect:/then: alone.
        // AgentRpcBroker.get() returns a process-wide singleton and a sibling spec builds one with
        // `tls: false`; if its cleanup ever regressed, this spec would inherit a cleartext broker and
        // fail with an opaque `Unsupported or unrecognized SSL message` out of startHandshake()
        // instead of naming the cause.
        assert registration.fingerprint ==~ /[0-9a-f]{64}/

        when: 'the certificate is taken off the socket, not out of the broker'
        def chain = servedChain(registration.endpoint)
        // computed here, not folded into the comparison, so a failure prints the two digests as the
        // bare 64-hex strings agent-rpc itself reports on a pinning rejection
        def presented = MessageDigest.getInstance('SHA-256').digest(chain[0].encoded).encodeHex().toString()

        then: 'the DER Netty wrote is byte for byte the DER the fingerprint was taken over'
        // The whole pinning design rests on this equality and nothing else asserts it: the task is
        // given a digest, never the certificate, so any re-encoding between
        // `builder.build(signer).getEncoded()` and the socket makes the two ends disagree forever.
        presented == registration.fingerprint

        and: 'the broker serves a single self-signed leaf'
        // Not because a chain would break pinning -- it would not, the leaf stays rawCerts[0] and the
        // digest still matches -- but because the broker deliberately serves one self-signed
        // certificate with no issuer above it. A chain appearing here means the credential shape
        // changed, and what the pin then commits to has to be re-examined rather than assumed.
        chain.length == 1
        ((X509Certificate) chain[0]).subjectX500Principal == ((X509Certificate) chain[0]).issuerX500Principal
    }

    /**
     * The certificate chain as RECEIVED from the broker over a real TLS socket.
     *
     * <p>The trust manager accepts anything, which is the proxy's own posture -- {@code agent-rpc}
     * sets {@code InsecureSkipVerify} and then decides by digest. Building a trust manager out of
     * {@code broker.certificatePem} instead would turn this into the same path-building check
     * {@code AgentRpcBrokerTest} already performs, i.e. exactly the comparison that cannot see a
     * re-encoding.
     *
     * <p>Deliberately NOT wrapped in a try/catch: a handshake that fails here IS the failure this
     * spec exists to report, and the only effect of catching it would be to keep the spec green
     * after the transport changed underneath it.
     */
    private static Certificate[] servedChain(String endpoint) {
        final parts = endpoint.split(':')
        final context = SSLContext.getInstance('TLS')
        context.init(null, [ new X509TrustManager() {
            @Override void checkClientTrusted(X509Certificate[] chain, String authType) {}
            @Override void checkServerTrusted(X509Certificate[] chain, String authType) {}
            @Override X509Certificate[] getAcceptedIssuers() { return new X509Certificate[0] }
        } ] as TrustManager[], null)
        final socket = (SSLSocket) context.socketFactory.createSocket(parts[0], parts[1] as int)
        try {
            socket.setSoTimeout(HANDSHAKE_TIMEOUT_MILLIS)
            // Offer h2 so the server runs the same ALPN selection it runs for the proxy. grpc's
            // server config is SelectedListenerFailureBehavior.ACCEPT, so an ALPN-less client would
            // also complete the handshake -- but then the negotiation production depends on would be
            // the one path this probe skipped, and the connection would be torn down by grpc's
            // "Failed protocol negotiation" branch instead of by the close below.
            final params = socket.getSSLParameters()
            params.setApplicationProtocols([ 'h2' ] as String[])
            socket.setSSLParameters(params)
            socket.startHandshake()
            return socket.session.peerCertificates
        }
        finally {
            // No h2 preface is ever sent, so the server sees a connection that completed TLS and then
            // went away. That is intentional: this spec is about the certificate bytes, and driving a
            // gRPC call here would only re-test what AgentRpcBrokerTest already covers.
            socket.close()
        }
    }
}
