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

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import software.amazon.awssdk.core.ResponseBytes
import software.amazon.awssdk.awscore.exception.AwsErrorDetails
import software.amazon.awssdk.awscore.exception.AwsServiceException
import software.amazon.awssdk.core.exception.SdkClientException
import software.amazon.awssdk.core.exception.SdkException
import nextflow.cloud.aws.AwsClientFactory
import software.amazon.awssdk.services.s3.model.CopyObjectRequest
import software.amazon.awssdk.services.s3.model.GetObjectRequest
import software.amazon.awssdk.services.s3.model.GetObjectResponse
import software.amazon.awssdk.services.s3.model.RequestPayer
import software.amazon.awssdk.services.s3.model.PutObjectRequest
import software.amazon.awssdk.services.s3.model.PutObjectResponse
import software.amazon.awssdk.services.s3.model.ServerSideEncryption
import software.amazon.awssdk.services.s3.model.NoSuchBucketException
import software.amazon.awssdk.services.s3.model.NoSuchKeyException
import software.amazon.awssdk.services.s3.model.Tag
import software.amazon.awssdk.services.s3.model.TaggingDirective
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for the AWS SDK → NIO exception conversion in {@link S3Client#convertAwsException}.
 */
class S3ClientTest extends Specification {

    def 'the claim PUT carries the bucket-policy settings its sibling writes apply'() {
        given: 'a client configured with server-side encryption, as an encrypting bucket requires'
        def sdk = Mock(software.amazon.awssdk.services.s3.S3Client)
        def factory = Mock(AwsClientFactory) { getS3Client(_, _) >> sdk }
        def client = new S3Client(factory, new Properties(), false)
        client.setStorageEncryption('aws:kms')
        client.setKmsKeyId('key-123')
        and:
        PutObjectRequest captured = null

        when:
        def created = client.putObjectIfAbsent('bkt', 'work/aa/bb/.command.claim')

        then:
        1 * sdk.putObject(_ as PutObjectRequest, _) >> { PutObjectRequest r, def body ->
            captured = r; return PutObjectResponse.builder().build()
        }
        created

        and: 'the conditional create is intact -- that is what makes the claim atomic'
        captured.ifNoneMatch() == '*'

        and: 'and SSE/KMS are applied, so a `deny unless encrypted` policy cannot 403 every claim'
        captured.serverSideEncryption() == ServerSideEncryption.AWS_KMS
        captured.ssekmsKeyId() == 'key-123'
    }

    def 'requester-pays reaches the claim PUT and the ranged GET, not only getObject'() {
        given: 'a client on a requester-pays bucket'
        def sdk = Mock(software.amazon.awssdk.services.s3.S3Client)
        def factory = Mock(AwsClientFactory) { getS3Client(_, _) >> sdk }
        def client = new S3Client(factory, new Properties(), false)
        client.setRequesterPaysEnabled('true')
        and:
        PutObjectRequest put = null
        GetObjectRequest get = null

        when: 'the work-dir claim -- without the payer this 403s, and a 403 is deliberately NOT read'
        and: 'as "lost the race", so the run dies with an error pointing nowhere near the cause'
        client.putObjectIfAbsent('bkt', 'work/aa/bb/.command.claim')
        then:
        1 * sdk.putObject(_ as PutObjectRequest, _) >> { PutObjectRequest r, def body ->
            put = r; return PutObjectResponse.builder().build()
        }
        put.requestPayer() == RequestPayer.REQUESTER

        when: 'the sampled identity\'s ranged read'
        client.getObjectRange(GetObjectRequest.builder().bucket('bkt').key('k').range('bytes=0-9').build())
        then:
        1 * sdk.getObjectAsBytes(_ as GetObjectRequest) >> { GetObjectRequest r ->
            get = r; return ResponseBytes.fromByteArray(GetObjectResponse.builder().build(), new byte[0])
        }
        get.requestPayer() == RequestPayer.REQUESTER
        and: 'the range the caller built is preserved'
        get.range() == 'bytes=0-9'

    }

    def 'should map NoSuchBucketException to NoSuchFileException'() {
        given:
        def aws = NoSuchBucketException.builder().message('nope').build()

        when:
        def result = S3Client.convertAwsException(aws, 'listObjects', 'my-bucket', null)

        then:
        result instanceof NoSuchFileException
        result.file == 's3://my-bucket'
        result.cause.is(aws)
    }

    def 'should map NoSuchKeyException to NoSuchFileException'() {
        given:
        def aws = NoSuchKeyException.builder().message('missing').build()

        when:
        def result = S3Client.convertAwsException(aws, 'getObject', 'my-bucket', 'path/to/obj')

        then:
        result instanceof NoSuchFileException
        result.file == 's3://my-bucket/path/to/obj'
        result.cause.is(aws)
    }

    @Unroll
    def 'should map HTTP #code to NoSuchFileException'() {
        given:
        def aws = AwsServiceException.builder()
                .message('err')
                .awsErrorDetails(AwsErrorDetails.builder().errorCode('X').build())
                .statusCode(code)
                .build()

        when:
        def result = S3Client.convertAwsException(aws, 'getObject', 'my-bucket', 'key')

        then:
        result instanceof NoSuchFileException
        result.file == 's3://my-bucket/key'
        result.cause.is(aws)

        where:
        code << [404]
    }

    @Unroll
    def 'should map HTTP #code to AccessDeniedException'() {
        given:
        def aws = AwsServiceException.builder()
                .message('denied')
                .awsErrorDetails(AwsErrorDetails.builder().errorCode('X').build())
                .statusCode(code)
                .build()

        when:
        def result = S3Client.convertAwsException(aws, 'getObject', 'my-bucket', 'key')

        then:
        result instanceof AccessDeniedException
        result.file == 's3://my-bucket/key'
        result.cause.is(aws)

        where:
        code << [401, 403]
    }

    def 'should map other AwsServiceException to generic IOException'() {
        given:
        def aws = AwsServiceException.builder()
                .message('boom')
                .awsErrorDetails(AwsErrorDetails.builder().errorCode('X').build())
                .statusCode(500)
                .build()

        when:
        def result = S3Client.convertAwsException(aws, 'putObject', 'my-bucket', 'k')

        then:
        result instanceof IOException
        !(result instanceof NoSuchFileException)
        !(result instanceof AccessDeniedException)
        result.message.contains('putObject')
        result.message.contains('s3://my-bucket/k')
        result.cause.is(aws)
    }

    def 'should map non-service SdkException to generic IOException'() {
        given:
        SdkException aws = SdkClientException.builder().message('network down').build()

        when:
        def result = S3Client.convertAwsException(aws, 'listBuckets', null, null)

        then:
        result.getClass() == IOException
        result.message.contains('listBuckets')
        result.message.contains('s3://')
        result.cause.is(aws)
    }

    @Unroll
    def 'should always set REPLACE tagging directive on copy so source tags are not inherited'() {
        given:
        def reqBuilder = CopyObjectRequest.builder()
                .sourceBucket('src').sourceKey('a')
                .destinationBucket('dst').destinationKey('b')

        when:
        S3Client.applyTagging(reqBuilder, tags)
        def req = reqBuilder.build()

        then:
        // REPLACE is always used so the destination never inherits the source object's tags
        req.taggingDirective() == TaggingDirective.REPLACE
        req.tagging() == expectedTagging

        where:
        tags                                            | expectedTagging
        null                                            | null
        []                                              | null
        [Tag.builder().key('foo').value('bar').build()] | 'foo=bar'
    }

    def 'should format path without trailing slash when key is null or empty'() {
        expect:
        (S3Client.convertAwsException(NoSuchBucketException.builder().message('').build(), 'op', 'b', key) as NoSuchFileException).file == 's3://b'

        where:
        key << [null, '']
    }
}
