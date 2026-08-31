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

import com.google.common.base.CaseFormat
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description

/**
 * Model the HTTP client timeout settings used to reach a remote service.
 * <p>
 * The two timeouts bound different phases: {@code connectTimeout} bounds establishing the
 * connection, {@code requestTimeout} bounds waiting for the answer once it is established.
 * Only the second one protects against a server that accepts the connection and then never
 * replies.
 * <p>
 * Values resolve the same way {@link RetryConfig} resolves its own: the config map first,
 * then the environment variable derived from {@code envPrefix} and the option name, then the
 * default. The env fallback matters because some clients — SCM repository access above all —
 * run before any Nextflow config has been parsed, so the config map is the only route that
 * exists for an embedder and the environment is the only route that exists for an operator.
 *
 * @author Rob Syme <rob.syme@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class HttpClientOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        The maximum time to wait for the connection to the remote service to be established,
        applied to a single attempt. It bounds only reaching the service, never waiting for
        its response — see `requestTimeout` for that. It must be greater than zero: unlike
        `requestTimeout`, `0 sec` is not accepted as "unbounded".
    """)
    final private Duration connectTimeout

    @ConfigOption
    @Description("""
        The maximum time to wait for a response from the remote service once the connection
        has been established, applied to a single attempt. Set to `0 sec` to wait
        indefinitely.
    """)
    final private Duration requestTimeout

    /**
     * @param config    The config map holding the {@code connectTimeout} / {@code requestTimeout}
     *                  entries, or an empty map to resolve from the environment
     * @param envPrefix The prefix of the environment variables to fall back to e.g. {@code NXF_GIT_}
     *                  yields {@code NXF_GIT_CONNECT_TIMEOUT} and {@code NXF_GIT_REQUEST_TIMEOUT}
     * @param defConnectTimeout The connect timeout to use when neither source provides a value
     * @param defRequestTimeout The request timeout to use when neither source provides a value
     */
    HttpClientOpts(Map config, String envPrefix, Duration defConnectTimeout, Duration defRequestTimeout) {
        connectTimeout = duration0(config, 'connectTimeout', envPrefix, defConnectTimeout)
        requestTimeout = duration0(config, 'requestTimeout', envPrefix, defRequestTimeout)
        // zero is not "unbounded" here — a connect phase that never gives up is not a bound
        // worth expressing, and HttpClient.Builder rejects it. Catch it at the config edge so
        // the message names the option rather than surfacing from the builder later on.
        if( connectTimeout.toMillis() <= 0 )
            throw new IllegalArgumentException("Config option 'connectTimeout' must be greater than zero - offending value: ${connectTimeout}")
    }

    /**
     * Resolve one option, reporting the source that carried the offending value when it does
     * not parse. {@link RetryConfig#valueOf} raises a bare "Not a valid duration value", which
     * is of little use when the value came from an environment variable the user has to find.
     */
    private static Duration duration0(Map config, String name, String envPrefix, Duration defValue) {
        try {
            return RetryConfig.valueOf(config, name, envPrefix, defValue, Duration)
        }
        catch( IllegalArgumentException e ) {
            final fromConfig = config?.get(name)
            final source = fromConfig != null
                ? "config option '${name}'"
                : "variable ${envKey(envPrefix, name)}"
            final value = fromConfig != null
                ? fromConfig
                : SysEnv.get(envKey(envPrefix, name))
            throw new IllegalArgumentException("Invalid duration value for ${source}: '${value}'", e)
        }
    }

    private static String envKey(String envPrefix, String name) {
        if( !envPrefix.endsWith('_') )
            envPrefix += '_'
        return envPrefix.toUpperCase() + CaseFormat.LOWER_CAMEL.to(CaseFormat.UPPER_UNDERSCORE, name)
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
        final long millis = requestTimeout.toMillis()
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
