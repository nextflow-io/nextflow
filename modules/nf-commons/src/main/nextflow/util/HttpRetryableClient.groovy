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
 *
 */

package nextflow.util

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.util.retry.Retryable

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.channels.UnresolvedAddressException
import java.time.Duration

/**
 * A utility class that provides HTTP client functionality with built-in retry capabilities
 * using the libseqera retry library.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class HttpRetryableClient {

    private static final List<Integer> HTTP_RETRYABLE_CODES = [429, 500, 502, 503, 504]
    
    private final HttpClient httpClient
    private final Retryable.Config retryConfig

    HttpRetryableClient(HttpClient httpClient, Retryable.Config retryConfig) {
        this.httpClient = httpClient
        this.retryConfig = retryConfig
    }

    /**
     * Sends an HTTP request with retry logic, using the provided response body handler
     */
    <T> HttpResponse<T> send(HttpRequest request, HttpResponse.BodyHandler<T> responseBodyHandler) {
        if (!retryConfig) {
            return httpClient.send(request, responseBodyHandler)
        }

        final retryable = Retryable.<HttpResponse<T>>of(retryConfig)
            .retryCondition((Throwable e) -> isRetryable(e))
            .retryIf((HttpResponse<T> response) -> isRetryableStatusCode(response))
            .onRetry({ event ->
                log.debug("HTTP request retry - attempt: ${event.attempt}, exception: ${event.failure?.message}, response: ${event.result}")
            })

        return retryable.apply(() -> httpClient.send(request, responseBodyHandler))
    }

    /**
     * Determines if an exception should trigger a retry
     */
    private boolean isRetryable(Throwable t) {
        // only retry SocketException and ignore generic IOException
        return t instanceof SocketException && !isCausedByUnresolvedAddressException(t)
    }

    /**
     * Determines if an HTTP response status code should trigger a retry
     */
    private boolean isRetryableStatusCode(HttpResponse<?> response) {
        return response.statusCode() in HTTP_RETRYABLE_CODES
    }

    /**
     * Checks if the exception is caused by an UnresolvedAddressException
     * which should not be retried as it's a permanent DNS failure
     */
    private boolean isCausedByUnresolvedAddressException(Throwable t) {
        if( t instanceof UnresolvedAddressException )
            return true
        if( t.cause==null )
            return false
        else
            return isCausedByUnresolvedAddressException(t.cause)
    }

    /**
     * Factory method to create an HttpRetryableClient with a given retry configuration
     */
    static HttpRetryableClient create(HttpClient httpClient, Retryable.Config retryConfig) {
        return new HttpRetryableClient(httpClient, retryConfig)
    }

    /**
     * Factory method to create an HttpRetryableClient with default configuration
     */
    static HttpRetryableClient create(HttpClient httpClient) {
        return new HttpRetryableClient(httpClient, null)
    }

    /**
     * Factory method to create an HttpRetryableClient with RetryConfig
     */
    static HttpRetryableClient create(HttpClient httpClient, RetryConfig retryConfig) {
        return new HttpRetryableClient(httpClient, retryConfig ? toRetryableConfig(retryConfig) : null)
    }

    /**
     * Converts RetryConfig to Retryable.Config
     */
    static Retryable.Config toRetryableConfig(RetryConfig config) {
        return new Retryable.Config() {
            Duration getDelay() { Duration.ofMillis(config.delay.toMillis()) }
            Duration getMaxDelay() { Duration.ofMillis(config.maxDelay.toMillis()) }
            int getMaxAttempts() { config.maxAttempts }
            double getJitter() { config.jitter }
            double getMultiplier() { 2.0d } // Default multiplier since RetryConfig doesn't have it
        }
    }
}