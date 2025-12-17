/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.cloud.aws.nio.util

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import software.amazon.awssdk.awscore.exception.AwsServiceException
import software.amazon.awssdk.core.exception.SdkClientException
import software.amazon.awssdk.services.s3.model.S3Exception
import spock.lang.Specification

/**
 * Tests for S3ObjectSummaryLookup exception translation
 *
 * @author Jonathan Manning
 */
class S3ObjectSummaryLookupTest extends Specification {

    S3ObjectSummaryLookup lookup = new S3ObjectSummaryLookup()

    def 'should translate 403 status to AccessDeniedException'() {
        given:
        def s3Exception = S3Exception.builder()
            .statusCode(403)
            .message("Access Denied")
            .build()

        when:
        def result = lookup.translateException(s3Exception, "s3://bucket/key")

        then:
        result instanceof AccessDeniedException
        result.message.contains("Access denied")
        result.message.contains("credentials")
    }

    def 'should translate 401 status to AccessDeniedException'() {
        given:
        def s3Exception = S3Exception.builder()
            .statusCode(401)
            .message("Unauthorized")
            .build()

        when:
        def result = lookup.translateException(s3Exception, "s3://bucket/key")

        then:
        result instanceof AccessDeniedException
        result.message.contains("Access denied")
    }

    def 'should translate 404 status to NoSuchFileException'() {
        given:
        def s3Exception = S3Exception.builder()
            .statusCode(404)
            .message("Not Found")
            .build()

        when:
        def result = lookup.translateException(s3Exception, "s3://bucket/key")

        then:
        result instanceof NoSuchFileException
    }

    def 'should translate credential loading errors to AccessDeniedException'() {
        given:
        def clientException = SdkClientException.builder()
            .message("Unable to load credentials from any of the providers in the chain")
            .build()

        when:
        def result = lookup.translateException(clientException, "s3://bucket/key")

        then:
        result instanceof AccessDeniedException
        result.message.contains("credentials")
    }

    def 'should translate marshall errors to AccessDeniedException'() {
        given:
        def clientException = SdkClientException.builder()
            .message("Unable to marshall request to JSON: Key cannot be empty")
            .build()

        when:
        def result = lookup.translateException(clientException, "s3://bucket/key")

        then:
        result instanceof AccessDeniedException
        result.message.contains("credentials")
    }

    def 'should wrap other SDK exceptions as IOException'() {
        given:
        def s3Exception = S3Exception.builder()
            .statusCode(500)
            .message("Internal Server Error")
            .build()

        when:
        def result = lookup.translateException(s3Exception, "s3://bucket/key")

        then:
        result instanceof IOException
        !(result instanceof AccessDeniedException)
        !(result instanceof NoSuchFileException)
        result.message.contains("s3://bucket/key")
        result.cause == s3Exception
    }
}
