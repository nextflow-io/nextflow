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

package nextflow.k8s.client

import java.security.KeyStore
import java.security.PrivateKey
import java.security.Signature
import java.security.cert.X509Certificate
import java.security.interfaces.ECPrivateKey

import org.bouncycastle.asn1.pkcs.PrivateKeyInfo
import org.bouncycastle.jce.provider.BouncyCastleProvider
import org.bouncycastle.openssl.PEMKeyPair
import org.bouncycastle.openssl.PEMParser
import org.bouncycastle.openssl.jcajce.JcaPEMKeyConverter
import spock.lang.Specification

class SSLUtilsTest extends Specification {

    static final String P256_PARAMS = '''\
        -----BEGIN EC PARAMETERS-----
        BggqhkjOPQMBBw==
        -----END EC PARAMETERS-----
        '''.stripIndent()

    static final String P256_SEC1 = '''\
        -----BEGIN EC PRIVATE KEY-----
        MHcCAQEEIEYFtFW0m33AUGc7uukjSzdoJHFr0P7ShcbIXVyY9MQEoAoGCCqGSM49
        AwEHoUQDQgAEASwIfdadalosP+hX6sLtbyzgOWztI/LdR13GY8lBVRVPJEE9an+d
        CvOir37f5sKhywuVAIM0s46siKjoDja8EA==
        -----END EC PRIVATE KEY-----
        '''.stripIndent()

    static final String P256_PKCS8 = '''\
        -----BEGIN PRIVATE KEY-----
        MIGHAgEAMBMGByqGSM49AgEGCCqGSM49AwEHBG0wawIBAQQgRgW0VbSbfcBQZzu6
        6SNLN2gkcWvQ/tKFxshdXJj0xAShRANCAAQBLAh91p1qWiw/6Ffqwu1vLOA5bO0j
        8t1HXcZjyUFVFU8kQT1qf50K86Kvft/mwqHLC5UAgzSzjqyIqOgONrwQ
        -----END PRIVATE KEY-----
        '''.stripIndent()

    static final String P256_CERT = '''\
        -----BEGIN CERTIFICATE-----
        MIIBczCCARmgAwIBAgIUcfB0NrKl+wTm7MSlXmzJb0vAcnkwCgYIKoZIzj0EAwIw
        DzENMAsGA1UEAwwEdGVzdDAeFw0yNjExMDIxOTI1MDVaFw0zNjEwMzAxOTI1MDVa
        MA8xDTALBgNVBAMMBHRlc3QwWTATBgcqhkjOPQIBBggqhkjOPQMBBwNCAAQBLAh9
        1p1qWiw/6Ffqwu1vLOA5bO0j8t1HXcZjyUFVFU8kQT1qf50K86Kvft/mwqHLC5UA
        gzSzjqyIqOgONrwQo1MwUTAdBgNVHQ4EFgQUugZZwpNhEszsv5wYlfAVrEvklpkw
        HwYDVR0jBBgwFoAUugZZwpNhEszsv5wYlfAVrEvklpkwDwYDVR0TAQH/BAUwAwEB
        /zAKBggqhkjOPQQDAgNIADBFAiEAk7zgVsUUPImn1Af8nXFipABObGbE1gGu6PdE
        ndmtzlYCIGid9AcNNIjWH+yT88aNTgniSPhA6QC+kkuRmASGOs9k
        -----END CERTIFICATE-----
        '''.stripIndent()

    static final String P384_SEC1 = '''\
        -----BEGIN EC PRIVATE KEY-----
        MIGkAgEBBDAIXZMO54qq/CH+jUY3ZOK5HVDbZuHvZCLme5mj3R3P3D2wo7TJXR6G
        q9jz3rTEFwSgBwYFK4EEACKhZANiAAS6xoZoQkHFpGnjJ9i8MKjbLTtdrzrowxT2
        sgpHsAN/saTIHSxhvi6X4vx51eHeBuEQwyht+VHWPKqyStP85Yv3CIPAYCvr3IU7
        MRaqta3oAe+m9TiFSp9EpvGIRHSJWFg=
        -----END EC PRIVATE KEY-----
        '''.stripIndent()

    static final String P384_CERT = '''\
        -----BEGIN CERTIFICATE-----
        MIIBsDCCATagAwIBAgIUT+qWUL7kIzj5Ei067kQmJVudSCQwCgYIKoZIzj0EAwIw
        DzENMAsGA1UEAwwEdGVzdDAeFw0yNjExMDIxOTI1MDVaFw0zNjEwMzAxOTI1MDVa
        MA8xDTALBgNVBAMMBHRlc3QwdjAQBgcqhkjOPQIBBgUrgQQAIgNiAAS6xoZoQkHF
        pGnjJ9i8MKjbLTtdrzrowxT2sgpHsAN/saTIHSxhvi6X4vx51eHeBuEQwyht+VHW
        PKqyStP85Yv3CIPAYCvr3IU7MRaqta3oAe+m9TiFSp9EpvGIRHSJWFijUzBRMB0G
        A1UdDgQWBBTSy67RzB5gdKS0zxieY1ePCeUiLDAfBgNVHSMEGDAWgBTSy67RzB5g
        dKS0zxieY1ePCeUiLDAPBgNVHRMBAf8EBTADAQH/MAoGCCqGSM49BAMCA2gAMGUC
        MDqfML+4izvl51Ib+La2GGvpcYZNPQMRnzOvx5VjwXKtJ7DqPgoC1B0DDg0XnC0m
        WAIxALJ/Tp93WQY28B1UN/EoN8QmUM/IsYaX0UonW4QQS0arFYnk26PJprETXV1F
        +ihFog==
        -----END CERTIFICATE-----
        '''.stripIndent()

    // SEC1 key without the optional [0] curve parameters
    static final String SEC1_NO_PARAMS = '''\
        -----BEGIN EC PRIVATE KEY-----
        MCUCAQEEIAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB
        -----END EC PRIVATE KEY-----
        '''.stripIndent()

    static final char[] PASS = 'secret'.toCharArray()

    private static KeyStore keyStore(String cert, String key) {
        SSLUtils.createKeyStore(
                new ByteArrayInputStream(cert.bytes),
                new ByteArrayInputStream(key.bytes),
                'EC', PASS, null, PASS)
    }

    private static boolean signAndVerify(KeyStore store) {
        final alias = store.aliases().nextElement()
        final key = (PrivateKey) store.getKey(alias, PASS)
        final cert = (X509Certificate) store.getCertificate(alias)
        final data = 'hello world'.bytes
        final signer = Signature.getInstance('SHA256withECDSA')
        signer.initSign(key)
        signer.update(data)
        final sig = signer.sign()
        final verifier = Signature.getInstance('SHA256withECDSA')
        verifier.initVerify(cert.publicKey)
        verifier.update(data)
        return verifier.verify(sig)
    }

    def 'should load EC key: #DESC' () {
        when:
        def store = keyStore(CERT, KEY)
        then:
        store.getKey(store.aliases().nextElement(), PASS) instanceof ECPrivateKey
        signAndVerify(store)

        where:
        DESC                    | CERT        | KEY
        'P-256 SEC1'            | P256_CERT   | P256_SEC1
        'P-256 SEC1 w/ params'  | P256_CERT   | P256_PARAMS + P256_SEC1
        'P-256 PKCS#8'          | P256_CERT   | P256_PKCS8
        'P-384 SEC1'            | P384_CERT   | P384_SEC1
    }

    def 'should fail on SEC1 key without curve parameters' () {
        when:
        keyStore(P256_CERT, SEC1_NO_PARAMS)
        then:
        def e = thrown(IOException)
        e.message.contains('curve parameters')
    }

    def 'should fail on unsupported EC key: #DESC' () {
        when:
        keyStore(P256_CERT, KEY)
        then:
        def e = thrown(IOException)
        e.message.startsWith('Unsupported EC algorithm')

        where:
        DESC                | KEY
        'wrong PEM type'    | '-----BEGIN RSA PRIVATE KEY-----\nAAAA\n-----END RSA PRIVATE KEY-----\n'
        'invalid PKCS#8'    | '-----BEGIN PRIVATE KEY-----\nMAMCAQA=\n-----END PRIVATE KEY-----\n'
    }

    /**
     * Read an EC key with Bouncy Castle, as SSLUtils did before it was removed.
     */
    private static ECPrivateKey bouncyCastleKey(String pem) {
        final converter = new JcaPEMKeyConverter().setProvider(new BouncyCastleProvider())
        final object = new PEMParser(new StringReader(pem)).readObject()
        if( object instanceof PEMKeyPair )
            return (ECPrivateKey) converter.getKeyPair(object).getPrivate()
        if( object instanceof PrivateKeyInfo )
            return (ECPrivateKey) converter.getPrivateKey(object)
        throw new IllegalArgumentException("Unexpected PEM object: ${object?.getClass()}")
    }

    def 'should read the same EC key as Bouncy Castle: #DESC' () {
        given:
        def store = keyStore(CERT, KEY)
        def key = (ECPrivateKey) store.getKey(store.aliases().nextElement(), PASS)
        def expected = bouncyCastleKey(REF)

        expect:
        key.s == expected.s
        key.params.curve == expected.params.curve
        key.params.generator == expected.params.generator
        key.params.order == expected.params.order
        key.params.cofactor == expected.params.cofactor
        // the JDK keeps the input PKCS#8 encoding while Bouncy Castle
        // re-encodes it, so the bytes only match for SEC1 keys
        !SAME_ENCODING || key.encoded == expected.encoded

        where:
        // Bouncy Castle could not read a leading EC PARAMETERS block,
        // so that case is compared against the bare SEC1 key
        DESC                    | CERT        | KEY                       | REF          | SAME_ENCODING
        'P-256 SEC1'            | P256_CERT   | P256_SEC1                 | P256_SEC1    | true
        'P-256 SEC1 w/ params'  | P256_CERT   | P256_PARAMS + P256_SEC1   | P256_SEC1    | true
        'P-256 PKCS#8'          | P256_CERT   | P256_PKCS8                | P256_PKCS8   | false
        'P-384 SEC1'            | P384_CERT   | P384_SEC1                 | P384_SEC1    | true
    }
}
