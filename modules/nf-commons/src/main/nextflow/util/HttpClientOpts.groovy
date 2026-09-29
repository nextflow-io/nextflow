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

package nextflow.util

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model the HTTP client timeout settings used to reach a remote service.
 * <p>
 * The two timeouts bound different phases: {@code connectTimeout} bounds establishing the
 * connection, {@code requestTimeout} bounds waiting for the answer once it is established.
 * Only the second one protects against a server that accepts the connection and then never
 * replies.
 * <p>
 * This is a plain value holder: it carries the two timeouts and derives the values the
 * {@link java.net.http.HttpClient} API expects. Resolving them from a config map or the
 * environment is left to the caller that owns the client — see e.g.
 * {@code nextflow.scm.RepositoryProvider} — so that each client keeps its own option names,
 * environment prefix and defaults, following the same resolution shape as {@link RetryConfig}.
 *
 * @author Rob Syme <rob.syme@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class HttpClientOpts {

    final private Duration connectTimeout
    final private Duration requestTimeout

    /**
     * @param connectTimeout The maximum time to wait for the connection to the remote service
     *                  to be established, applied to a single attempt. It must be greater than
     *                  zero: unlike {@code requestTimeout}, {@code 0} is not accepted as "unbounded".
     * @param requestTimeout The maximum time to wait for a response once the connection has been
     *                  established, applied to a single attempt. {@code 0} means wait indefinitely.
     */
    HttpClientOpts(Duration connectTimeout, Duration requestTimeout) {
        // zero is not "unbounded" here — a connect phase that never gives up is not a bound
        // worth expressing, and HttpClient.Builder rejects it anyway
        if( connectTimeout==null || connectTimeout.toMillis() <= 0 )
            throw new IllegalArgumentException("Config option 'connectTimeout' must be greater than zero - offending value: ${connectTimeout}")
        this.connectTimeout = connectTimeout
        this.requestTimeout = requestTimeout
    }

    /**
     * The per-attempt connect timeout, in the form {@link java.net.http.HttpClient} expects.
     *
     * @return The timeout; never null, and guaranteed strictly positive
     */
    java.time.Duration connectTimeout() {
        return java.time.Duration.ofMillis(connectTimeout.toMillis())
    }

    /**
     * The per-attempt response timeout.
     *
     * @return The timeout, or null to wait indefinitely
     */
    java.time.Duration requestTimeout() {
        final long millis = requestTimeout != null ? requestTimeout.toMillis() : 0
        return millis > 0 ? java.time.Duration.ofMillis(millis) : null
    }

    /**
     * The value to pass to {@link java.net.http.HttpRequest.Builder#timeout}.
     * <p>
     * That timeout spans the connect phase as well as the response, so handing it the bare
     * {@link #requestTimeout} would cap the connection attempt at a value that may be shorter
     * than {@link #connectTimeout} — and a connect stall would then surface as
     * {@code HttpTimeoutException} rather than {@code HttpConnectTimeoutException}, i.e. as
     * the fail-fast case rather than the retryable one. Summing the two keeps each phase
     * bounded by its own setting.
     *
     * @return The timeout, or null when the request timeout is unbounded
     */
    java.time.Duration httpRequestTimeout() {
        final request = requestTimeout()
        return request != null
            ? connectTimeout().plus(request)
            : null
    }
}
