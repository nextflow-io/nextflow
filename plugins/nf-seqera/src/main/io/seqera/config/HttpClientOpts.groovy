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

package io.seqera.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model the HTTP client settings used to reach the Seqera scheduler service.
 * <p>
 * The two timeouts bound different phases: {@code connectTimeout} bounds reaching the
 * scheduler, {@code requestTimeout} bounds waiting for its answer. Only the second one
 * protects against a stalled scheduler, which is the case that hangs the task monitor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class HttpClientOpts implements ConfigScope {

    static final Duration DEFAULT_CONNECT_TIMEOUT = Duration.of('10 sec')

    static final Duration DEFAULT_REQUEST_TIMEOUT = Duration.of('45 sec')

    @ConfigOption
    @Description("""
        The maximum time to wait for a connection to the Seqera scheduler service to be
        established, applied to a single attempt (default: `10 sec`). It bounds only reaching
        the service, never waiting for its response — see `requestTimeout` for that. A
        timed-out connection proves the request never arrived, so it is retried under
        `seqera.executor.retryPolicy` for every request, task submissions included. It must be
        greater than zero — unlike `requestTimeout`, `0 sec` is not accepted as "unbounded".
    """)
    final private Duration connectTimeout

    @ConfigOption
    @Description("""
        The maximum time to wait for a response from the Seqera scheduler service, applied to
        a single attempt (default: `45 sec`). A timed-out poll or cancel is retried under
        `seqera.executor.retryPolicy`, so the worst case for one call is `requestTimeout`
        times `retryPolicy.maxAttempts` plus the accumulated delays; a timed-out task
        submission is not re-sent, since the scheduler may already be creating the tasks.
        Set to `0 sec` to wait indefinitely.
    """)
    final private Duration requestTimeout

    HttpClientOpts() {
        this(Collections.emptyMap())
    }

    HttpClientOpts(Map opts) {
        // an explicit null check rather than Groovy truthiness: Duration.asBoolean() is
        // false at zero, so `?:` would quietly turn the documented "0 sec means unbounded"
        // back into the default
        requestTimeout = opts.requestTimeout != null
            ? opts.requestTimeout as Duration
            : DEFAULT_REQUEST_TIMEOUT
        connectTimeout = opts.connectTimeout != null
            ? opts.connectTimeout as Duration
            : DEFAULT_CONNECT_TIMEOUT
        // unlike requestTimeout, zero is not "unbounded" here — a connect phase that never
        // gives up is not a bound worth expressing, and sched-client rejects it. Catch it at
        // the config edge so the message names the option rather than surfacing from the
        // client builder halfway through session start.
        if( connectTimeout.toMillis() <= 0 )
            throw new IllegalArgumentException("Config option 'seqera.executor.httpClient.connectTimeout' must be greater than zero - offending value: ${opts.connectTimeout}")
    }

    /**
     * The per-attempt connect timeout, in the form {@code sched-client} expects.
     *
     * @return the timeout; never null, and rejected upstream unless strictly positive
     */
    java.time.Duration connectTimeout() {
        return java.time.Duration.ofMillis(connectTimeout.toMillis())
    }

    /**
     * The per-attempt response timeout, in the form {@code sched-client} expects.
     * <p>
     * The timeout bounds a single attempt, and what it bounds overall depends on the request.
     * A stalled poll or cancel is idempotent, so the client re-sends it: the timeout raises an
     * {@code IOException} the retry policy absorbs, and the task is retried rather than failed
     * out of {@code checkAllTasks()}. The cost is that the retries run on the
     * {@code Task monitor} thread, so a scheduler stalled on every attempt holds that thread
     * for up to {@code requestTimeout × retryPolicy.maxAttempts} plus the accumulated backoff
     * — bounded, where before it was not, but far from free. A stalled task submission is not
     * re-sent, so it is bounded at a single timeout however many attempts remain.
     *
     * @return the timeout, or null to wait indefinitely — the {@code sched-client} default
     */
    java.time.Duration requestTimeout() {
        final long millis = requestTimeout.toMillis()
        return millis > 0 ? java.time.Duration.ofMillis(millis) : null
    }
}
