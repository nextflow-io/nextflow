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

import software.amazon.awssdk.awscore.exception.AwsErrorDetails
import software.amazon.awssdk.awscore.exception.AwsServiceException
import software.amazon.awssdk.core.exception.SdkClientException
import software.amazon.awssdk.core.exception.SdkException
import software.amazon.awssdk.services.s3.model.NoSuchBucketException
import software.amazon.awssdk.services.s3.model.NoSuchKeyException
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for the AWS SDK → NIO exception conversion in {@link S3Client#convertAwsException}.
 */
class S3ClientTest extends Specification {

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

    def 'should format path without trailing slash when key is null or empty'() {
        expect:
        (S3Client.convertAwsException(NoSuchBucketException.builder().message('').build(), 'op', 'b', key) as NoSuchFileException).file == 's3://b'

        where:
        key << [null, '']
    }
}
