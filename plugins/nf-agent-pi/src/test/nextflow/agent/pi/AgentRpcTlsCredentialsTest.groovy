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
import java.security.Security
import java.security.cert.CertificateFactory
import java.security.cert.X509Certificate

import spock.lang.Specification

/**
 * Guards the two properties of the broker's TLS identity that fail silently or globally if broken:
 * the fingerprint is the digest of the certificate DER (the value the Go proxy compares), and
 * generating it registers no JVM-wide JCA provider.
 */
class AgentRpcTlsCredentialsTest extends Specification {

    private static X509Certificate parse(String pem) {
        return (X509Certificate) CertificateFactory.getInstance('X.509')
            .generateCertificate(new ByteArrayInputStream(pem.bytes))
    }

    def 'should fingerprint the certificate DER as lowercase hex'() {
        when:
        def credentials = AgentRpcTlsCredentials.create()
        def certificate = parse(credentials.certificatePem)

        then: 'the digest is over getEncoded() -- the DER -- not over the PEM text or its base64 body'
        credentials.fingerprint ==~ /[0-9a-f]{64}/
        credentials.fingerprint == MessageDigest.getInstance('SHA-256')
            .digest(certificate.encoded).encodeHex().toString()

        and: 'the PEM the digest was taken from is what a TLS stack will parse back'
        credentials.certificateStream().text == credentials.certificatePem
        credentials.privateKeyStream().text.startsWith('-----BEGIN PRIVATE KEY-----')
    }

    def 'should generate a self-signed EC certificate usable as a server identity'() {
        when:
        def credentials = AgentRpcTlsCredentials.create()
        def certificate = parse(credentials.certificatePem)

        then: 'EC P-256 and self-signed: the pin makes the issuer irrelevant'
        certificate.publicKey.algorithm == 'EC'
        certificate.sigAlgName == 'SHA256withECDSA'
        certificate.subjectX500Principal == certificate.issuerX500Principal
        certificate.subjectX500Principal.name.contains('nextflow-agent-rpc')

        and: 'valid now, with room for clock skew and a long pipeline'
        certificate.checkValidity()

        and: 'loopback names, so a conventional TLS client can also verify it'
        certificate.subjectAlternativeNames*.last() as Set == ['localhost', '127.0.0.1'] as Set
    }

    def 'should mint a distinct identity per broker and touch no global JCA state'() {
        given:
        def before = Security.providers*.name

        when:
        def first = AgentRpcTlsCredentials.create()
        def second = AgentRpcTlsCredentials.create()

        then: 'a fresh key pair each time -- the identity is per run, never reused across runs'
        first.fingerprint != second.fingerprint
        first.certificatePem != second.certificatePem

        and: 'no Security.addProvider: a second plugin must not mutate the JVM provider list'
        Security.providers*.name == before
        Security.getProvider('BC') == null
    }
}
