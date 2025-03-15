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

package io.seqera.client

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration
import java.time.temporal.ChronoUnit
import java.util.concurrent.Executors
import java.util.function.Predicate

import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.config.SeqeraConfig
import io.seqera.exception.BadResponseException
import io.seqera.sched.api.v1a1.CancelJobResponse
import io.seqera.sched.api.v1a1.CreateJobRequest
import io.seqera.sched.api.v1a1.CreateJobResponse
import io.seqera.sched.api.v1a1.DescribeJobResponse
import io.seqera.sched.api.v1a1.GetJobLogsResponse
import io.seqera.serde.Serde
import io.seqera.util.trace.TraceUtils
import nextflow.Session
import nextflow.util.Threads
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SeqeraClient {

    private HttpClient httpClient

    private CookieManager cookieManager

    private SeqeraConfig config

    SeqeraClient(Session session) {
        this.config = new SeqeraConfig(session.config.seqera as Map ?: Collections.emptyMap())
        // the cookie manager
        cookieManager = new CookieManager()
        // create http client
        this.httpClient = newHttpClient()
    }

    protected HttpClient newHttpClient() {
        final builder = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_1_1)
            .followRedirects(HttpClient.Redirect.NORMAL)
            .connectTimeout(Duration.ofSeconds(10))
            .cookieHandler(cookieManager)
        // use virtual threads executor if enabled
        if( Threads.useVirtual() )
            builder.executor(Executors.newVirtualThreadPerTaskExecutor())
        // build and return the new client
        return builder.build()
    }

    protected <R> R get0(String uri, Class<R> responseType) {
        final trace = TraceUtils.rndTrace()
        final req = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .headers('Content-Type','application/json', 'Traceparent', trace)
            .GET()
            .build()
        send0(req, responseType)
    }

    protected <R> R delete0(String uri, Class<R> response) {
        final trace = TraceUtils.rndTrace()
        final req = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .headers('Content-Type','application/json', 'Traceparent', trace)
            .DELETE()
            .build()
        return send0(req, response)
    }

    protected <T,R> R post0(String uri, T request, Class<R> response) {
        final trace = TraceUtils.rndTrace()
        final body = JsonOutput.toJson(request)
        final req = HttpRequest.newBuilder()
            .uri(URI.create(uri))
            .headers('Content-Type','application/json', 'Traceparent', trace)
            .POST(HttpRequest.BodyPublishers.ofString(body))
            .build()
        return send0(req, response)
    }

    protected <T> T send0(HttpRequest req, Class<T> type) {
        final HttpResponse<String> resp = httpSend(req);
        log.debug("Seqera status response: statusCode={}; body={}", resp.statusCode(), resp.body())
        if( resp.statusCode()==200 ) {
            return Serde.fromJson(resp.body(), type)
        }
        else {
            final msg = "Seqera invalid response: ${req.method()} ${req.uri()} [${resp.statusCode()}] ${resp.body()}"
            throw new BadResponseException(msg)
        }
    }

    CreateJobResponse createJob(CreateJobRequest request) {
        final uri = "${config.endpoint}/v1a1/jobs"
        return post0(uri, request, CreateJobResponse)
    }

    DescribeJobResponse describeJob(String jobId) {
        final uri = "${config.endpoint}/v1a1/jobs/${jobId}"
        return get0(uri, DescribeJobResponse)
    }

    CancelJobResponse cancelJob(String jobId) {
        final uri = "${config.endpoint}/v1a1/jobs/${jobId}"
        return delete0(uri, CancelJobResponse)
    }

    GetJobLogsResponse getJobLogs(String jobId) {
        final uri = "${config.endpoint}/v1a1/jobs/${jobId}/logs"
        return get0(uri, GetJobLogsResponse)
    }

    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond, Predicate<T> handle) {
        final cfg = config.retryOpts()
        final listener = new EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                def msg = "Seqera connection failure - attempt: ${event.attemptCount}"
                if( event.lastResult!=null )
                    msg += "; response: ${event.lastResult}"
                if( event.lastFailure != null )
                    msg += "; exception: [${event.lastFailure.class.name}] ${event.lastFailure.message}"
                log.debug(msg)
            }
        }
        return RetryPolicy.<T>builder()
            .handleIf(cond)
            .handleResultIf(handle)
            .withBackoff(cfg.delay.toMillis(), cfg.maxDelay.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(cfg.maxAttempts)
            .withJitter(cfg.jitter)
            .onRetry(listener)
            .build()
    }

    protected <T> HttpResponse<T> safeApply(CheckedSupplier action) {
        final retryOnException = (e -> e instanceof IOException) as Predicate<? extends Throwable>
        final retryOnStatusCode = ((HttpResponse<T> resp) -> resp.statusCode() in SERVER_ERRORS) as Predicate<HttpResponse<T>>
        final policy = retryPolicy(retryOnException, retryOnStatusCode)
        return Failsafe.with(policy).get(action)
    }

    static private final List<Integer> SERVER_ERRORS = [429,500,502,503,504]

    protected HttpResponse<String> httpSend(HttpRequest req)  {
        return safeApply(() -> httpClient.send(req, HttpResponse.BodyHandlers.ofString()))
    }
}
