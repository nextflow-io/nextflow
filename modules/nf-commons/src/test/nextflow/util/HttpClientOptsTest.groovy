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

import spock.lang.Specification

/**
 *
 * @author Rob Syme <rob.syme@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HttpClientOptsTest extends Specification {

    def 'should expose the connect and request timeouts as java.time durations' () {
        given:
        def opts = new HttpClientOpts(Duration.of('10s'), Duration.of('250ms'))

        expect:
        opts.connectTimeout() == java.time.Duration.ofSeconds(10)
        opts.requestTimeout() == java.time.Duration.ofMillis(250)
    }

    def 'should treat a zero or missing request timeout as unbounded' () {
        expect:
        new HttpClientOpts(Duration.of('10s'), Duration.of('0s')).requestTimeout() == null
        new HttpClientOpts(Duration.of('10s'), null).requestTimeout() == null
    }

    def 'should reject a non-positive connect timeout' () {
        when:
        new HttpClientOpts(connect, Duration.of('30s'))
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('connectTimeout')
        e.message.contains('greater than zero')

        where:
        connect << [Duration.of('0s'), null]
    }

    def 'should sum connect and request for the http request timeout, or none when unbounded' () {
        expect:
        // HttpRequest.timeout() spans the connect phase as well as the response, so a
        // request timeout shorter than the connect timeout would surface a connect stall
        // as the non-retryable request-timeout exception
        new HttpClientOpts(Duration.of('10s'), Duration.of('45s')).httpRequestTimeout() == java.time.Duration.ofSeconds(55)
        and: 'an unbounded request timeout applies none'
        new HttpClientOpts(Duration.of('10s'), Duration.of('0s')).httpRequestTimeout() == null
    }
}
