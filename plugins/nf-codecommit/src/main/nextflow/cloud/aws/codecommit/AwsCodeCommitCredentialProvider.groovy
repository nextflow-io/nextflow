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
 */

package nextflow.cloud.aws.codecommit

import javax.crypto.Mac

/*
 * Copyright 2013-2019 the original author or authors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import javax.crypto.spec.SecretKeySpec
import javax.security.auth.login.CredentialException
import java.security.MessageDigest
import java.security.NoSuchAlgorithmException
import java.text.SimpleDateFormat

import com.amazonaws.auth.AWSCredentials
import com.amazonaws.auth.AWSCredentialsProvider
import com.amazonaws.auth.AWSSessionCredentials
import com.amazonaws.auth.AWSStaticCredentialsProvider
import com.amazonaws.auth.BasicAWSCredentials
import com.amazonaws.auth.DefaultAWSCredentialsProviderChain
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import org.eclipse.jgit.errors.UnsupportedCredentialItem
import org.eclipse.jgit.transport.CredentialItem
import org.eclipse.jgit.transport.CredentialsProvider
import org.eclipse.jgit.transport.URIish
/**
 * Provides a jgit {@link CredentialsProvider} implementation that can provide the
 * appropriate credentials to connect to an AWS CodeCommit repository.
 *
 * From the command line, you can configure git to use AWS CodeCommit with a credential
 * helper. Although jgit does not support credential helper commands, it provides
 * a CredentialsProvider abstract class we can extend. Connecting to an AWS CodeCommit
 * (codecommit) repository requires an AWS access key and secret key. These are used to
 * calculate a signature for the git request. The AWS access key is used as the codecommit
 * username, and the calculated signature is used as the password. The process for
 * calculating this signature is documented at
 * https://docs.aws.amazon.com/general/latest/gr/signature-version-4.html.
 *
 * @author Don Laidlaw
 *
 * This class is mostly a direct port of the AwsCodeCommitCredentialProvider class provided by the
 * spring framework in org.springframwork.cloud.config.server.support, with some minor
 * modifications and simplifications that leverage the Groovy language.
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 */
@Slf4j
@CompileStatic
final class AwsCodeCommitCredentialProvider extends CredentialsProvider {

    private static final String SHA_256 = "SHA-256"
    private static final String UTF8 = "UTF8"
    private static final String HMAC_SHA256 = "HmacSHA256"

    private String username
    private String password

    /**
     * Calculate the AWS CodeCommit password for the provided URI and AWS secret key. This
     * uses the algorithm published by AWS at
     * https://docs.aws.amazon.com/general/latest/gr/signature-version-4.html
     * @param uri the codecommit repository uri
     * @param awsSecretKey the aws secret key
     * @return the password to use in the git request
     */
    protected static String calculateCodeCommitPassword(URIish uri, String awsSecretKey, Date now) {
        String[] split = uri.getHost().split("\\.")
        if (split.length < 4) {
            throw new CredentialException("Cannot detect AWS region from URI $uri")
        }
        String region = split[1]

        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd'T'HHmmss")
        dateFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

        String dateStamp = dateFormat.format(now)
        String shortDateStamp = dateStamp.substring(0, 8)

        String codeCommitPassword
        try {
            def stringToSign = "AWS4-HMAC-SHA256\n${dateStamp}\n${shortDateStamp}/$region/codecommit/aws4_request\n${bytesToHexString(canonicalRequestDigest(uri))}"
            def signedRequest = bytesToHexString(
                    sign(awsSecretKey, shortDateStamp, region, stringToSign)
            )

            codeCommitPassword = "${dateStamp}Z${signedRequest}"
        }
        catch (Exception e) {
            throw new RuntimeException("Error calculating AWS CodeCommit password", e)
        }

        return codeCommitPassword
    }

    private static byte[] hmacSha256(String data, byte[] key) {
        String algorithm = HMAC_SHA256
        Mac mac = Mac.getInstance(algorithm)
        mac.init(new SecretKeySpec(key, algorithm))
        return mac.doFinal(data.getBytes(UTF8))
    }

    private static byte[] sign(String secret, String shortDateStamp, String region, String toSign) {
        byte[] kSecret = ("AWS4" + secret).getBytes(UTF8)
        byte[] kDate = hmacSha256(shortDateStamp, kSecret)
        byte[] kRegion = hmacSha256(region, kDate)
        byte[] kService = hmacSha256("codecommit", kRegion)
        byte[] kSigning = hmacSha256("aws4_request", kService)
        return hmacSha256(toSign, kSigning)
    }

    /**
     * Creates a message digest.
     * @param uri uri to process
     * @return a message digest
     * @throws NoSuchAlgorithmException when the SHA 256 algorithm is not found
     */
    private static byte[] canonicalRequestDigest(URIish uri) {
        def canonicalRequest = ""

        // this could be done faster with a templated multi-line GString, but is broken out here
        // to maintain documentation of each part
        canonicalRequest += "GIT\n"  						// codecommit uses GIT as the request method
        canonicalRequest += "${uri.getPath()}\n"  			// URI request path
        canonicalRequest += "\n"  							// Query string, always empty for codecommit

        // Next are canonical headers — codecommit only requires the host header
        canonicalRequest += "host:${uri.getHost()}\n\n"  	// canonical headers are alays terminated with \n
        canonicalRequest += "host\n"  						// The list of canonical headers, only one for codecommit

        MessageDigest digest = MessageDigest.getInstance(SHA_256);

        return digest.digest(canonicalRequest.getBytes());
    }

    /**
     * Convert bytes to a hex string.
     * @param bytes the bytes
     * @return a string of hex characters encoding the bytes.
     */
    private static String bytesToHexString(byte[] bytes) { bytes.encodeHex().toString() }

    /**
     * This provider can handle uris like
     * https://git-codecommit.$AWS_REGION.amazonaws.com/v1/repos/$REPO .
     * @param uri uri to parse
     * @return {@code true} if the URI can be handled
     */
    static boolean canHandle(String uri) {
        if ( !uri?.trim() ) {
            return false
        }

        try {
            URL url = new URL(uri)
            URI u = new URI(url.getProtocol(), url.getUserInfo(), url.getHost(),
                    url.getPort(), url.getPath(), url.getQuery(), url.getRef())
            if (u.getScheme().equals("https")) {
                String host = u.getHost()
                if (host.endsWith(".amazonaws.com")
                        && host.startsWith("git-codecommit.")) {
                    return true
                }
            }
        }
        catch (Throwable t) {
            // ignore all, we can't handle it
            log.debug "AWS CodeCommit cannot handle uri: $uri - Reason: ${t.message ?: t}"
        }

        return false
    }

    /**
     * Get the AWSCredentials. If an AWSCredentialProvider was specified, use that,
     * otherwise, create a new AWSCredentialsProvider. If the username and password are
     * provided, use those directly as AWSCredentials. Otherwise, use the
     * {@link DefaultAWSCredentialsProviderChain} as is standard with AWS applications.
     * @return the AWS credentials.
     */
    private AWSCredentials retrieveAwsCredentials() {
        AWSCredentialsProvider credsProvider
        if ( username && password ) {
            log.debug "Creating a static AWS credentials provider"
            credsProvider = new AWSStaticCredentialsProvider( new BasicAWSCredentials( username, password ))
        }
        else {
            log.debug "Creating a default AWS credentials provider chain"
            credsProvider = new DefaultAWSCredentialsProviderChain()
        }
        return credsProvider.getCredentials()
    }

    /**
     * This credentials provider cannot run interactively.
     * @return false
     * @see org.eclipse.jgit.transport.CredentialsProvider#isInteractive()
     */
    @Override
    boolean isInteractive() { false }

    /**
     * We support username and password credential items only.
     * @see org.eclipse.jgit.transport.CredentialsProvider#supports(org.eclipse.jgit.transport.CredentialItem[])
     */
    @Override
    boolean supports(CredentialItem... items) {
        for ( i in items ) {
            if (i instanceof CredentialItem.Username) {
                continue
            }
            else if (i instanceof CredentialItem.Password) {
                continue
            }
            else {
                return false
            }
        }
        return true
    }

    /**
     * Get the username and password to use for the given uri.
     * @see org.eclipse.jgit.transport.CredentialsProvider#get(org.eclipse.jgit.transport.URIish,
     * org.eclipse.jgit.transport.CredentialItem[])
     */
    @Override
    boolean get(URIish uri, CredentialItem... items) {
        String codeCommitPassword
        String awsAccessKey
        String awsSecretKey

        try {
            AWSCredentials awsCredentials = retrieveAwsCredentials()
            StringBuilder awsKey = new StringBuilder();
            awsKey.append(awsCredentials.getAWSAccessKeyId());
            awsSecretKey = awsCredentials.getAWSSecretKey();
            if (awsCredentials instanceof AWSSessionCredentials) {
                AWSSessionCredentials sessionCreds = (AWSSessionCredentials) awsCredentials;
                if ( sessionCreds.getSessionToken() ) {
                    awsKey.append('%').append(sessionCreds.getSessionToken())
                }
            }
            awsAccessKey = awsKey.toString()
        }
        catch (Exception e) {
            throw new AbortOperationException("Unable to retrieve AWS Credentials", e)
        }

        try {
            codeCommitPassword = calculateCodeCommitPassword(uri, awsSecretKey, new Date());
        }
        catch (Exception e) {
            throw new AbortOperationException("Error calculating AWS CodeCommit password", e)
        }

        for ( i in items ) {
            if (i instanceof CredentialItem.Username) {
                ((CredentialItem.Username) i).setValue(awsAccessKey);
                log.trace("Returning username " + awsAccessKey);
                continue;
            }
            if (i instanceof CredentialItem.Password) {
                ((CredentialItem.Password) i).setValue(codeCommitPassword.toCharArray());
                log.trace("Returning password " + codeCommitPassword);
                continue;
            }
            if (i instanceof CredentialItem.StringType
                    && i.getPromptText().equals("Password: ")) {
                ((CredentialItem.StringType) i).setValue(codeCommitPassword);
                log.trace("Returning password string " + codeCommitPassword);
                continue;
            }
            throw new UnsupportedCredentialItem(uri, i.getClass().getName() + ":" + i.getPromptText());
        }

        return true;
    }

    /**
     * Throw out cached data and force retrieval of AWS credentials.
     * @param uri This parameter is not used in this implementation.
     */
    @Override
    void reset(URIish uri) {
        // Should throw out cached info.
        // Note that even though the credentials (password) we calculate here is
        // valid for 15 minutes, we do not cache it. Instead, we re-calculate
        // it each time it is needed. However, the AWSCredentialProvider will cache
        // its AWSCredentials object.
    }

    /**
     * @param awsCredentialProvider the awsCredentialProvider to set
     */
    void setAwsCredentialsProvider(AWSCredentialsProvider awsCredentialsProvider) {
        this.awsCredentialsProvider = awsCredentialsProvider
    }

    String getUsername() {
        return username
    }

    String getPassword() {
        return password
    }

}
