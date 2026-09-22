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

package nextflow.cloud.aws.nio

import static com.github.tomakehurst.wiremock.client.WireMock.*

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.core.WireMockConfiguration
import nextflow.cloud.aws.AwsClientFactory
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.model.GetObjectRequest
import spock.lang.Specification

/**
 * The claim PUT asks the SDK to retry a 409 {@code ConditionalRequestConflict}, which the SDK does
 * not retry by default. That retry happens INSIDE the SDK client, so a mocked client cannot reach
 * it -- these tests stub the HTTP layer underneath instead and count the requests that arrive.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class S3ClientRetryTest extends Specification {

    WireMockServer wireMock

    def setup() {
        wireMock = new WireMockServer(WireMockConfiguration.options().dynamicPort())
        wireMock.start()
        // the S3Client constructor probes the caller account; answer it so it does not add noise
        wireMock.stubFor(get(urlMatching('/\\?.*')).willReturn(aResponse().withStatus(404)))
    }

    def cleanup() {
        wireMock?.stop()
    }

    /** An {@code S3Client} whose underlying SDK client really talks HTTP, to the stub above. */
    private S3Client clientOf() {
        final sdk = software.amazon.awssdk.services.s3.S3Client.builder()
                .endpointOverride(URI.create(wireMock.baseUrl()))
                .region(Region.US_EAST_1)
                .credentialsProvider(StaticCredentialsProvider.create(AwsBasicCredentials.create('k', 's')))
                .forcePathStyle(true)
                .build()
        return new S3Client(Mock(AwsClientFactory) { getS3Client(_, _) >> sdk }, new Properties(), false)
    }

    private static String errorXml(String code) {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?><Error><Code>${code}</Code>" +
                "<Message>${code}</Message><RequestId>req-1</RequestId></Error>"
    }

    def 'a 409 conflict is retried, and the claim is won when the retry succeeds'() {
        given: 'the first conditional PUT conflicts with a concurrent one, the retry gets through'
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('conflict').whenScenarioStateIs('Started')
                .willReturn(aResponse().withStatus(409).withBody(errorXml('ConditionalRequestConflict')))
                .willSetStateTo('retried'))
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('conflict').whenScenarioStateIs('retried')
                .willReturn(aResponse().withStatus(200).withHeader('ETag', '"abc"')))

        when:
        def created = clientOf().putObjectIfAbsent('bkt', 'lock')

        then: 'this caller created the object -- which mapping the 409 to `false` would have denied'
        created
        and: 'the SDK really did re-issue the request'
        wireMock.verify(2, putRequestedFor(urlEqualTo('/bkt/lock')))
    }

    def 'a 409 conflict is retried, and the claim is lost when the retry finds the object'() {
        given: 'the competing PUT landed in between, so the retry sees the object already there'
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('conflict').whenScenarioStateIs('Started')
                .willReturn(aResponse().withStatus(409).withBody(errorXml('ConditionalRequestConflict')))
                .willSetStateTo('retried'))
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('conflict').whenScenarioStateIs('retried')
                .willReturn(aResponse().withStatus(412).withBody(errorXml('PreconditionFailed'))))

        when:
        def created = clientOf().putObjectIfAbsent('bkt', 'lock')

        then: 'the retry resolved the ambiguity the 409 left: the other writer won'
        !created
        wireMock.verify(2, putRequestedFor(urlEqualTo('/bkt/lock')))
    }

    def 'OperationAborted is a 409 too, and is retried on the same footing'() {
        given: 'the other 409 S3 answers a conflicting conditional operation with'
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('aborted').whenScenarioStateIs('Started')
                .willReturn(aResponse().withStatus(409).withBody(errorXml('OperationAborted')))
                .willSetStateTo('retried'))
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .inScenario('aborted').whenScenarioStateIs('retried')
                .willReturn(aResponse().withStatus(200).withHeader('ETag', '"abc"')))

        when:
        def created = clientOf().putObjectIfAbsent('bkt', 'lock')

        then: 'matching the STATUS covers it -- an errorCode check for ConditionalRequestConflict'
        and: 'alone would have missed this one and reported a claim that was never attempted'
        created
        wireMock.verify(2, putRequestedFor(urlEqualTo('/bkt/lock')))
    }

    def 'a 412 is definitive and is not retried'() {
        given:
        wireMock.stubFor(put(urlEqualTo('/bkt/lock'))
                .willReturn(aResponse().withStatus(412).withBody(errorXml('PreconditionFailed'))))

        when:
        def created = clientOf().putObjectIfAbsent('bkt', 'lock')

        then: 'the object exists -- there is nothing to resolve, so no second request'
        !created
        wireMock.verify(1, putRequestedFor(urlEqualTo('/bkt/lock')))
    }

    def 'without the per-request plugin the SDK does not retry a 409 at all'() {
        given: 'the same status on a call that does NOT add the retry predicate'
        wireMock.stubFor(get(urlEqualTo('/bkt/obj'))
                .willReturn(aResponse().withStatus(409).withBody(errorXml('ConditionalRequestConflict'))))

        when:
        clientOf().getObjectRange(GetObjectRequest.builder().bucket('bkt').key('obj').range('bytes=0-9').build())

        then: 'it throws on the first response -- this is what makes the plugin above load-bearing'
        thrown(IOException)
        wireMock.verify(1, getRequestedFor(urlEqualTo('/bkt/obj')))
    }
}
