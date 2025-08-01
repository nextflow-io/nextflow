/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.util

import io.seqera.util.retry.Retryable
import spock.lang.Specification
import spock.lang.Unroll

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.channels.UnresolvedAddressException
import java.time.Duration

/**
 * Unit tests for HttpRetryableClient
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class HttpRetryableClientTest extends Specification {

    def "should create client with retry config"() {
        given:
        def httpClient = Mock(HttpClient)
        def retryConfig = Mock(Retryable.Config)

        when:
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        then:
        client != null
        client.httpClient == httpClient
        client.retryConfig == retryConfig
    }

    def "should create client without retry config"() {
        given:
        def httpClient = Mock(HttpClient)

        when:
        def client = HttpRetryableClient.create(httpClient)

        then:
        client != null
        client.httpClient == httpClient
        client.retryConfig == null
    }

    def "should create client with RetryConfig"() {
        given:
        def httpClient = Mock(HttpClient)
        def retryConfig = Mock(RetryConfig) {
            getDelay() >> nextflow.util.Duration.of('500ms')
            getMaxDelay() >> nextflow.util.Duration.of('30s')
            getMaxAttempts() >> 5
            getJitter() >> 0.25d
        }

        when:
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        then:
        client != null
        client.httpClient == httpClient
        client.retryConfig != null
        client.retryConfig.delay == Duration.ofMillis(500)
        client.retryConfig.maxDelay == Duration.ofMillis(30000)
        client.retryConfig.maxAttempts == 5
        client.retryConfig.jitter == 0.25d
        client.retryConfig.multiplier == 2.0d
    }

    def "should send request without retry when no config"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def bodyHandler = HttpResponse.BodyHandlers.ofString()
        def response = Mock(HttpResponse)
        def client = HttpRetryableClient.create(httpClient)

        when:
        def result = client.send(request, bodyHandler)

        then:
        1 * httpClient.send(request, bodyHandler) >> response
        result == response
    }

    def "should send request successfully"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def response = Mock(HttpResponse) {
            statusCode() >> 200
            body() >> "success"
        }
        def retryConfig = Mock(RetryConfig) {
            getDelay() >> nextflow.util.Duration.of('500ms')
            getMaxDelay() >> nextflow.util.Duration.of('30s')
            getMaxAttempts() >> 5
            getJitter() >> 0.25d
        }

        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        def result = client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> response
        0 * httpClient.send(request, HttpResponse.BodyHandlers.ofString())
        result == response
    }

    @Unroll
    def "should handle retriable HTTP status code: #statusCode"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 2
            getDelay() >> Duration.ofMillis(10)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)
        def response = Mock(HttpResponse)

        when:
        def result = client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> response
        1 * response.statusCode() >> statusCode
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> response
        1 * response.statusCode() >> statusCode
        0 * httpClient.send(request, HttpResponse.BodyHandlers.ofString())
        result == response

        where:
        statusCode << [429, 500, 502, 503, 504]
    }

    def "should handle SocketException with retries"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 2  // Only 1 attempt to avoid complex mocking
            getDelay() >> Duration.ofMillis(10)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        2 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> {
            throw new SocketException("Connection refused") 
        }
        thrown(SocketException)
    }

    def "should not retry on UnresolvedAddressException"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 3
            getDelay() >> Duration.ofMillis(100)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> { 
            throw new UnresolvedAddressException() 
        }
        thrown(UnresolvedAddressException)
    }

    def "should not retry on SocketException caused by UnresolvedAddressException"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 3
            getDelay() >> Duration.ofMillis(100)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)
        def unresolvedCause = new UnresolvedAddressException()
        def socketException = new SocketException("Connection failed")
        socketException.initCause(unresolvedCause)

        when:
        client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> { throw socketException }
        thrown(SocketException)
    }

    @Unroll
    def "should not retry on non-retryable HTTP status code: #statusCode"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 3
            getDelay() >> Duration.ofMillis(100)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def response = Mock(HttpResponse) {
            statusCode() >> statusCode
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        def result = client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> response
        result == response

        where:
        statusCode << [200, 201, 301, 400, 401, 403, 404]
    }

    def "should not retry on IOException that is not SocketException"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 3
            getDelay() >> Duration.ofMillis(100)
            getMaxDelay() >> Duration.ofSeconds(1)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        1 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> { 
            throw new IOException("Generic IO error") 
        }
        thrown(IOException)
    }


    def "should convert RetryConfig to Retryable.Config correctly"() {
        given:
        def retryConfig = Mock(RetryConfig) {
            getDelay() >> nextflow.util.Duration.of('1s')
            getMaxDelay() >> nextflow.util.Duration.of('2m')
            getMaxAttempts() >> 10
            getJitter() >> 0.5d
        }

        when:
        def retryableConfig = HttpRetryableClient.toRetryableConfig(retryConfig)

        then:
        retryableConfig.delay == Duration.ofSeconds(1)
        retryableConfig.maxDelay == Duration.ofMinutes(2)
        retryableConfig.maxAttempts == 10
        retryableConfig.jitter == 0.5d
        retryableConfig.multiplier == 2.0d
    }

    def "should exhaust retries and throw exception"() {
        given:
        def httpClient = Mock(HttpClient)
        def request = Mock(HttpRequest)
        def retryConfig = Mock(Retryable.Config) {
            getMaxAttempts() >> 2
            getDelay() >> Duration.ofMillis(10)
            getMaxDelay() >> Duration.ofMillis(100)
            getJitter() >> 0.1d
            getMultiplier() >> 2.0d
        }
        def client = HttpRetryableClient.create(httpClient, retryConfig)

        when:
        client.send(request, HttpResponse.BodyHandlers.ofString())

        then:
        2 * httpClient.send(request, HttpResponse.BodyHandlers.ofString()) >> { 
            throw new SocketException("Always fails") 
        }
        thrown(SocketException)
    }
}