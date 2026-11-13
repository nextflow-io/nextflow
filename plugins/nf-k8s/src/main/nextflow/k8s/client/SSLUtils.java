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

package nextflow.k8s.client;

/**
 * This file is derived from
 * https://github.com/kubernetes-client/java/blob/master/util/src/main/java/io/kubernetes/client/util/SSLUtils.java
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.security.KeyFactory;
import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.PrivateKey;
import java.security.UnrecoverableKeyException;
import java.security.cert.Certificate;
import java.security.cert.CertificateException;
import java.security.cert.CertificateFactory;
import java.security.cert.X509Certificate;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.PKCS8EncodedKeySpec;
import java.security.spec.RSAPrivateCrtKeySpec;
import javax.net.ssl.KeyManager;
import javax.net.ssl.KeyManagerFactory;

import org.apache.commons.codec.binary.Base64;


public class SSLUtils {

    // DER encoding of OID 1.2.840.10045.2.1 (id-ecPublicKey)
    private static final byte[] EC_PUBLIC_KEY_OID = {0x06, 0x07, 0x2A, (byte) 0x86, 0x48, (byte) 0xCE, 0x3D, 0x02, 0x01};

    public static boolean isNotNullOrEmpty(String val) {
        return val != null && val.length() > 0;
    }

    public static KeyManager[] keyManagers(String certData, String certFile, String keyData, String keyFile,
                                           String algo, String passphrase, String keyStoreFile, String keyStorePassphrase)
            throws NoSuchAlgorithmException, UnrecoverableKeyException, KeyStoreException, CertificateException,
            InvalidKeySpecException, IOException {
        KeyManager[] keyManagers = null;
        if ((isNotNullOrEmpty(certData) || isNotNullOrEmpty(certFile))
                && (isNotNullOrEmpty(keyData) || isNotNullOrEmpty(keyFile))) {
            KeyStore keyStore = createKeyStore(certData, certFile, keyData, keyFile, algo, passphrase, keyStoreFile,
                    keyStorePassphrase);
            KeyManagerFactory kmf = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm());
            kmf.init(keyStore, passphrase.toCharArray());
            keyManagers = kmf.getKeyManagers();
        }
        return keyManagers;
    }


    public static KeyStore createKeyStore(String clientCertData, String clientCertFile, String clientKeyData,
                                          String clientKeyFile, String clientKeyAlgo, String clientKeyPassphrase, String keyStoreFile,
                                          String keyStorePassphrase) throws IOException, CertificateException, NoSuchAlgorithmException,
            InvalidKeySpecException, KeyStoreException {
        try (InputStream certInputStream = getInputStreamFromDataOrFile(clientCertData, clientCertFile);
             InputStream keyInputStream = getInputStreamFromDataOrFile(clientKeyData, clientKeyFile)) {
            return createKeyStore(certInputStream, keyInputStream, clientKeyAlgo,
                    clientKeyPassphrase != null ? clientKeyPassphrase.toCharArray() : null,
                    keyStoreFile, getKeyStorePassphrase(keyStorePassphrase));
        }
    }

    /**
     * Read an EC private key in PKCS#8 (BEGIN PRIVATE KEY) or SEC1 (BEGIN EC PRIVATE KEY) PEM format.
     */
    static private PrivateKey generateEcKey(InputStream keyInputStream) throws IOException {
        PemBlock pem = readPem(keyInputStream);
        byte[] pkcs8;
        if( "PRIVATE KEY".equals(pem.type) )
            pkcs8 = pem.bytes;
        else if( "EC PRIVATE KEY".equals(pem.type) )
            pkcs8 = sec1ToPkcs8(pem.bytes);
        else
            throw new IOException("Unsupported EC algorithm: unexpected PEM block '" + pem.type + "'");
        try {
            return KeyFactory.getInstance("EC").generatePrivate(new PKCS8EncodedKeySpec(pkcs8));
        }
        catch (NoSuchAlgorithmException | InvalidKeySpecException e) {
            throw new IOException("Unsupported EC algorithm: " + e.getMessage(), e);
        }
    }

    /**
     * Wrap a SEC1 ECPrivateKey into a PKCS#8 PrivateKeyInfo, taking the
     * named curve from the SEC1 [0] parameters field.
     */
    private static byte[] sec1ToPkcs8(byte[] sec1) throws IOException {
        Asn1Object sequence = new DerParser(sec1).read();
        sequence.validateSequence();
        DerParser parser = new DerParser(sequence.getValue());
        if( !BigInteger.ONE.equals(parser.read().getInteger()) )
            throw new IOException("Unsupported EC algorithm: invalid SEC1 key version");
        parser.read();
        byte[] curveOid = null;
        while( parser.hasMore() ) {
            Asn1Object obj = parser.read();
            if( obj.tag == 0xA0 ) {
                curveOid = obj.getValue();
                break;
            }
        }
        if( curveOid == null )
            throw new IOException("Unsupported EC algorithm: SEC1 key does not specify the curve parameters");
        if( curveOid.length == 0 || curveOid[0] != 0x06 )
            throw new IOException("Unsupported EC algorithm: only named curves are supported");

        byte[] algorithm = derEncode(0x30, EC_PUBLIC_KEY_OID, curveOid);
        byte[] version = {0x02, 0x01, 0x00};
        return derEncode(0x30, version, algorithm, derEncode(0x04, sec1));
    }

    private static byte[] derEncode(int tag, byte[]... parts) {
        ByteArrayOutputStream content = new ByteArrayOutputStream();
        for( byte[] part : parts )
            content.writeBytes(part);
        int length = content.size();
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        out.write(tag);
        if( length < 0x80 ) {
            out.write(length);
        }
        else {
            byte[] lenBytes = BigInteger.valueOf(length).toByteArray();
            int offset = lenBytes[0] == 0 ? 1 : 0;
            out.write(0x80 | (lenBytes.length - offset));
            out.write(lenBytes, offset, lenBytes.length - offset);
        }
        out.writeBytes(content.toByteArray());
        return out.toByteArray();
    }

    static private PrivateKey generateStdKey(InputStream keyInputStream, String clientKeyAlgo) throws IOException, NoSuchAlgorithmException, InvalidKeySpecException {
        byte[] keyBytes = decodePem(keyInputStream);

        KeyFactory keyFactory = KeyFactory.getInstance(clientKeyAlgo);
        try {
            // First let's try PKCS8
            return keyFactory.generatePrivate(new PKCS8EncodedKeySpec(keyBytes));
        }
        catch (InvalidKeySpecException e) {
            // Otherwise try PKCS1
            RSAPrivateCrtKeySpec keySpec = decodePKCS1(keyBytes);
            return keyFactory.generatePrivate(keySpec);
        }
    }

    public static KeyStore createKeyStore(InputStream certInputStream, InputStream keyInputStream, String clientKeyAlgo,
                                          char[] clientKeyPassphrase, String keyStoreFile, char[] keyStorePassphrase) throws IOException,
            CertificateException, NoSuchAlgorithmException, InvalidKeySpecException, KeyStoreException {
        CertificateFactory certFactory = CertificateFactory.getInstance("X509");
        X509Certificate cert = (X509Certificate) certFactory.generateCertificate(certInputStream);

        PrivateKey privateKey  = clientKeyAlgo.equals("EC")
                    ? generateEcKey(keyInputStream)
                    : generateStdKey(keyInputStream, clientKeyAlgo);

        KeyStore keyStore = KeyStore.getInstance("JKS");
        if (keyStoreFile != null && keyStoreFile.length() > 0) {
            keyStore.load(new FileInputStream(keyStoreFile), keyStorePassphrase);
        } else {
            loadDefaultKeyStoreFile(keyStore, keyStorePassphrase);
        }

        String alias = cert.getSubjectX500Principal().getName();
        keyStore.setKeyEntry(alias, privateKey, clientKeyPassphrase, new Certificate[] { cert });

        return keyStore;
    }

    static class PemBlock {
        final String type;
        final byte[] bytes;

        PemBlock(String type, byte[] bytes) {
            this.type = type;
            this.bytes = bytes;
        }
    }

    private static byte[] decodePem(InputStream keyInputStream) throws IOException {
        return readPem(keyInputStream).bytes;
    }

    // This method is inspired and partly taken over from
    // http://oauth.googlecode.com/svn/code/java/
    // All credits to belong to them.
    private static PemBlock readPem(InputStream keyInputStream) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(keyInputStream));
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.contains("-----BEGIN ")) {
                    String begin = line.trim();
                    String type = begin.replace("-----BEGIN ", "").replace("-----", "");
                    byte[] bytes = readBytes(reader, begin.replace("BEGIN", "END"));
                    // openssl ecparam -genkey emits the curve parameters before the key
                    if( "EC PARAMETERS".equals(type) )
                        continue;
                    return new PemBlock(type, bytes);
                }
            }
            throw new IOException("PEM is invalid: no begin marker");
        } finally {
            reader.close();
        }
    }

    private static byte[] readBytes(BufferedReader reader, String endMarker) throws IOException {
        String line;
        StringBuffer buf = new StringBuffer();

        while ((line = reader.readLine()) != null) {
            if (line.indexOf(endMarker) != -1) {
                return Base64.decodeBase64(buf.toString());
            }
            buf.append(line.trim());
        }
        throw new IOException("PEM is invalid : No end marker");
    }

    public static RSAPrivateCrtKeySpec decodePKCS1(byte[] keyBytes) throws IOException {
        DerParser parser = new DerParser(keyBytes);
        Asn1Object sequence = parser.read();
        sequence.validateSequence();
        parser = new DerParser(sequence.getValue());
        parser.read();

        return new RSAPrivateCrtKeySpec(next(parser), next(parser), next(parser), next(parser), next(parser),
                next(parser), next(parser), next(parser));
    }

    private static BigInteger next(DerParser parser) throws IOException {
        return parser.read().getInteger();
    }

    static class DerParser {

        private InputStream in;

        DerParser(byte[] bytes) throws IOException {
            this.in = new ByteArrayInputStream(bytes);
        }

        Asn1Object read() throws IOException {
            int tag = in.read();

            if (tag == -1) {
                throw new IOException("Invalid DER: stream too short, missing tag");
            }

            int length = getLength();
            byte[] value = new byte[length];
            if (in.read(value) < length) {
                throw new IOException("Invalid DER: stream too short, missing value");
            }

            return new Asn1Object(tag, value);
        }

        boolean hasMore() throws IOException {
            return in.available() > 0;
        }

        private int getLength() throws IOException {
            int i = in.read();
            if (i == -1) {
                throw new IOException("Invalid DER: length missing");
            }

            if ((i & ~0x7F) == 0) {
                return i;
            }

            int num = i & 0x7F;
            if (i >= 0xFF || num > 4) {
                throw new IOException("Invalid DER: length field too big (" + i + ")");
            }

            byte[] bytes = new byte[num];
            if (in.read(bytes) < num) {
                throw new IOException("Invalid DER: length too short");
            }

            return new BigInteger(1, bytes).intValue();
        }
    }

    static class Asn1Object {

        private final int type;
        private final byte[] value;
        private final int tag;

        public Asn1Object(int tag, byte[] value) {
            this.tag = tag;
            this.type = tag & 0x1F;
            this.value = value;
        }

        public byte[] getValue() {
            return value;
        }

        BigInteger getInteger() throws IOException {
            if (type != 0x02) {
                throw new IOException("Invalid DER: object is not integer"); //$NON-NLS-1$
            }
            return new BigInteger(value);
        }

        void validateSequence() throws IOException {
            if (type != 0x10) {
                throw new IOException("Invalid DER: not a sequence");
            }
            if ((tag & 0x20) != 0x20) {
                throw new IOException("Invalid DER: can't parse primitive entity");
            }
        }
    }

    private static void loadDefaultKeyStoreFile(KeyStore keyStore, char[] keyStorePassphrase)
            throws CertificateException, NoSuchAlgorithmException, IOException {

        String keyStorePath = System.getProperty("javax.net.ssl.keyStore");
        if (keyStorePath != null && keyStorePath.length() > 0) {
            File keyStoreFile = new File(keyStorePath);
            if (loadDefaultStoreFile(keyStore, keyStoreFile, keyStorePassphrase)) {
                return;
            }
        }

        keyStore.load(null);
    }

    private static boolean loadDefaultStoreFile(KeyStore keyStore, File fileToLoad, char[] passphrase)
            throws CertificateException, NoSuchAlgorithmException, IOException {
        if (fileToLoad.exists() && fileToLoad.isFile() && fileToLoad.length() > 0) {
            keyStore.load(new FileInputStream(fileToLoad), passphrase);
            return true;
        }
        return false;
    }

    public static InputStream getInputStreamFromDataOrFile(String data, String file) throws FileNotFoundException {
        if (data != null) {
            byte[] bytes = Base64.decodeBase64(data);
            // TODO handle non-base64 here?
            return new ByteArrayInputStream(bytes);
        }
        if (file != null) {
            return new FileInputStream(file);
        }
        return null;
    }

    private static char[] getKeyStorePassphrase(String keyStorePassphrase) {
        if (keyStorePassphrase == null || keyStorePassphrase.length() == 0) {
            return System.getProperty("javax.net.ssl.keyStorePassword", "changeit").toCharArray();
        }
        return keyStorePassphrase.toCharArray();
    }
}
