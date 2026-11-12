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

import java.nio.charset.StandardCharsets
import java.security.KeyPair
import java.security.KeyPairGenerator
import java.security.MessageDigest
import java.security.SecureRandom
import java.security.spec.ECGenParameterSpec
import java.util.concurrent.TimeUnit

import groovy.transform.CompileStatic
import org.bouncycastle.asn1.x500.X500Name
import org.bouncycastle.asn1.x509.Extension
import org.bouncycastle.asn1.x509.GeneralName
import org.bouncycastle.asn1.x509.GeneralNames
import org.bouncycastle.cert.jcajce.JcaX509v3CertificateBuilder
import org.bouncycastle.operator.jcajce.JcaContentSignerBuilder

/**
 * The driver-side TLS identity of the agent RPC broker: a key pair and a matching self-signed
 * certificate, generated once per run and held only in memory.
 *
 * <p>There is no CA, no trust store and no certificate file. An agent task is given the
 * certificate's SHA-256 {@link #getFingerprint fingerprint} on its command line and pins it — the
 * SSH known-hosts pattern — so the certificate needs no issuer anyone else trusts, and its subject
 * name is never what authenticates the driver. A fingerprint is a public commitment, not a secret,
 * so it is harmless in {@code .command.sh}.
 *
 * <p>Two properties are load-bearing and must not be broken casually:
 *
 * <ul>
 * <li><b>The fingerprint is the SHA-256 of the certificate's DER encoding, lowercase hex, 64
 *     characters, no separators and no algorithm prefix.</b> {@code agent-rpc}'s
 *     {@code VerifyPeerCertificate} hashes {@code rawCerts[0]}, which is exactly that DER. Hashing
 *     the PEM text instead — or its base64 body — yields an equally well-formed digest that simply
 *     never matches, and the resulting failure is indistinguishable from a genuine pinning
 *     rejection. Both ends must agree byte for byte.</li>
 * <li><b>No global JCA provider is registered.</b> Certificate construction needs only
 *     BouncyCastle's builder classes; the signature and the key pair come from the default
 *     provider. {@code nf-k8s} does call {@code Security.addProvider}, and a second plugin
 *     mutating JVM-wide state is exactly what this avoids.</li>
 * </ul>
 *
 * <p>EC P-256 rather than RSA: every agent task holds an open stream for its whole life, so the
 * per-connection handshake cost is what matters, not the one-off generation.
 */
@CompileStatic
class AgentRpcTlsCredentials {

    private static final String SUBJECT = 'CN=nextflow-agent-rpc'

    /** Signed with the key it certifies -- pinning makes the issuer irrelevant. */
    private static final String SIGNATURE_ALGORITHM = 'SHA256withECDSA'

    /**
     * Certificate validity. The identity dies with the run either way, so this only has to be wide
     * enough that a long pipeline and a modest clock skew between driver and execution node cannot
     * expire it mid-flight.
     */
    private static final long BACKDATE_MILLIS = TimeUnit.HOURS.toMillis(1)
    private static final long VALIDITY_MILLIS = TimeUnit.DAYS.toMillis(365)

    private final String certificatePem
    private final String privateKeyPem
    private final String fingerprint

    private AgentRpcTlsCredentials(String certificatePem, String privateKeyPem, String fingerprint) {
        this.certificatePem = certificatePem
        this.privateKeyPem = privateKeyPem
        this.fingerprint = fingerprint
    }

    /** Generates a fresh key pair and self-signed certificate. */
    static AgentRpcTlsCredentials create() {
        final random = new SecureRandom()
        final generator = KeyPairGenerator.getInstance('EC')
        generator.initialize(new ECGenParameterSpec('secp256r1'), random)
        final KeyPair keys = generator.generateKeyPair()
        final byte[] certificate = certificate(keys, random)
        return new AgentRpcTlsCredentials(
            pem('CERTIFICATE', certificate),
            pem('PRIVATE KEY', keys.private.encoded),
            fingerprint(certificate))
    }

    private static byte[] certificate(KeyPair keys, SecureRandom random) {
        final subject = new X500Name(SUBJECT)
        final now = System.currentTimeMillis()
        final builder = new JcaX509v3CertificateBuilder(
            subject,
            new BigInteger(64, random),
            new Date(now - BACKDATE_MILLIS),
            new Date(now + VALIDITY_MILLIS),
            subject,
            keys.public)
        // Deliberately no basicConstraints: a certificate that asserts CA=false cannot serve as its
        // own trust anchor for Go's verifier, which would break a client that later chose to pin by
        // trust root instead of by digest. Absent is permissive; CA=false is a trap.
        //
        // The names below are ones a pinning client never checks, but they keep the certificate
        // usable with a conventional TLS client (openssl s_client, this broker's tests) on loopback.
        builder.addExtension(Extension.subjectAlternativeName, false, new GeneralNames([
            new GeneralName(GeneralName.dNSName, 'localhost'),
            new GeneralName(GeneralName.iPAddress, '127.0.0.1') ] as GeneralName[]))
        final signer = new JcaContentSignerBuilder(SIGNATURE_ALGORITHM).build(keys.private)
        return builder.build(signer).getEncoded()
    }

    /**
     * SHA-256 over the certificate DER, lowercase hex. This is the exact string the proxy compares
     * against the digest of the certificate it is served -- see the class note.
     */
    private static String fingerprint(byte[] der) {
        return MessageDigest.getInstance('SHA-256').digest(der).encodeHex().toString()
    }

    private static String pem(String type, byte[] der) {
        final body = Base64.getMimeEncoder(64, '\n'.getBytes(StandardCharsets.US_ASCII)).encodeToString(der)
        return "-----BEGIN ${type}-----\n${body}\n-----END ${type}-----\n".toString()
    }

    String getCertificatePem() { certificatePem }

    String getFingerprint() { fingerprint }

    InputStream certificateStream() { new ByteArrayInputStream(certificatePem.getBytes(StandardCharsets.US_ASCII)) }

    InputStream privateKeyStream() { new ByteArrayInputStream(privateKeyPem.getBytes(StandardCharsets.US_ASCII)) }
}
