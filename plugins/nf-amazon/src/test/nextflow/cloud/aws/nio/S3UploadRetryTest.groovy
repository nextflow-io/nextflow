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

import java.net.URI
import java.nio.ByteBuffer

import nextflow.cloud.aws.nio.util.ByteBufferInputStream
import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.core.sync.RequestBody
import software.amazon.awssdk.http.AbortableInputStream
import software.amazon.awssdk.http.ExecutableHttpRequest
import software.amazon.awssdk.http.HttpExecuteRequest
import software.amazon.awssdk.http.HttpExecuteResponse
import software.amazon.awssdk.http.SdkHttpClient
import software.amazon.awssdk.http.SdkHttpResponse
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.PutObjectRequest
import spock.lang.Specification

/**
 * A single-part {@code S3OutputStream} upload sends its buffer as
 * {@code RequestBody.fromInputStream(new ByteBufferInputStream(buf), len)}. The SDK asks that body
 * for a fresh stream on EVERY attempt, so the upload is only retryable if the stream can be re-read.
 *
 * This drives a real {@link S3Client} over a stubbed transport that fails the first attempt with
 * 503, which is the shape of the transient error the SDK is meant to absorb.
 */
class S3UploadRetryTest extends Specification {

    /** Transport stub: reads the request body like a real send, fails attempt 1 with 503. */
    static class FlakyTransport implements SdkHttpClient {
        int attempts = 0
        List<Integer> bytesSeen = []

        @Override
        ExecutableHttpRequest prepareRequest(HttpExecuteRequest request) {
            return new ExecutableHttpRequest() {
                @Override
                HttpExecuteResponse call() {
                    attempts++
                    // consume the body, as sending it would
                    final provider = request.contentStreamProvider()
                    bytesSeen << (provider.present ? provider.get().newStream().bytes.length : 0)
                    final status = attempts == 1 ? 503 : 200
                    return HttpExecuteResponse.builder()
                            .response(SdkHttpResponse.builder().statusCode(status).build())
                            .responseBody(AbortableInputStream.create(new ByteArrayInputStream(new byte[0])))
                            .build()
                }

                @Override
                void abort() { }
            }
        }

        @Override
        void close() { }
    }

    def 'a single-part upload survives a transient 503 -- the body is re-read on the retry'() {
        given:
        def payload = ('x' * 128).bytes
        def transport = new FlakyTransport()
        def s3 = S3Client.builder()
                .region(Region.US_EAST_1)
                .credentialsProvider(StaticCredentialsProvider.create(AwsBasicCredentials.create('key', 'secret')))
                .endpointOverride(URI.create('https://s3.local'))
                .httpClient(transport)
                .build()
        and: 'the body S3OutputStream builds for a single-part putObject'
        def body = RequestBody.fromInputStream(
                new ByteBufferInputStream(ByteBuffer.wrap(payload)), payload.length)

        when:
        s3.putObject(PutObjectRequest.builder().bucket('bucket').key('key').build(), body)

        then: 'the SDK retried, and the retry re-sent as many bytes as the first attempt'
        transport.attempts == 2
        transport.bytesSeen[1] == transport.bytesSeen[0]
        and: 'each attempt carried the whole payload (plus the SDK aws-chunked framing)'
        transport.bytesSeen[0] >= payload.length

        cleanup:
        s3?.close()
    }
}
