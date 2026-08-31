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

import nextflow.SysEnv
import spock.lang.Specification

/**
 *
 * @author Rob Syme <rob.syme@gmail.com>
 */
class HttpClientOptsTest extends Specification {

    static final Duration CONNECT_DEF = Duration.of('60s')
    static final Duration REQUEST_DEF = Duration.of('45s')

    private HttpClientOpts opts(Map config=[:]) {
        new HttpClientOpts(config, 'NXF_GIT_', CONNECT_DEF, REQUEST_DEF)
    }

    def 'should use the given defaults when nothing is configured' () {
        given:
        SysEnv.push([:])

        expect:
        opts().connectTimeout() == java.time.Duration.ofSeconds(60)
        opts().requestTimeout() == java.time.Duration.ofSeconds(45)

        cleanup:
        SysEnv.pop()
    }

    def 'should take the value from the config map' () {
        given:
        SysEnv.push([:])

        expect:
        opts([connectTimeout: '5s']).connectTimeout() == java.time.Duration.ofSeconds(5)
        opts([requestTimeout: '250ms']).requestTimeout() == java.time.Duration.ofMillis(250)

        cleanup:
        SysEnv.pop()
    }

    def 'should fall back to the env variable when the config map has no value' () {
        given:
        SysEnv.push([NXF_GIT_CONNECT_TIMEOUT: '5s', NXF_GIT_REQUEST_TIMEOUT: '250ms'])

        expect:
        opts().connectTimeout() == java.time.Duration.ofSeconds(5)
        opts().requestTimeout() == java.time.Duration.ofMillis(250)

        cleanup:
        SysEnv.pop()
    }

    def 'should prefer the config map over the env variable' () {
        given:
        SysEnv.push([NXF_GIT_CONNECT_TIMEOUT: '5s'])

        expect:
        opts([connectTimeout: '30s']).connectTimeout() == java.time.Duration.ofSeconds(30)

        cleanup:
        SysEnv.pop()
    }

    def 'should treat a zero request timeout as unbounded' () {
        given:
        SysEnv.push([:])

        expect:
        opts([requestTimeout: '0s']).requestTimeout() == null

        cleanup:
        SysEnv.pop()
    }

    def 'should reject a non-positive connect timeout' () {
        given:
        SysEnv.push([:])

        when:
        opts([connectTimeout: '0s'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('connectTimeout')
        e.message.contains('greater than zero')

        cleanup:
        SysEnv.pop()
    }

    def 'should report the offending env variable when the value is not a duration' () {
        given:
        SysEnv.push([NXF_GIT_REQUEST_TIMEOUT: 'foo'])

        when:
        opts()
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('NXF_GIT_REQUEST_TIMEOUT')
        e.message.contains('foo')

        cleanup:
        SysEnv.pop()
    }

    def 'should sum connect and request for the http request timeout' () {
        given:
        SysEnv.push([:])

        expect:
        // HttpRequest.timeout() spans the connect phase as well as the response, so a
        // request timeout shorter than the connect timeout would surface a connect stall
        // as the non-retryable read-timeout exception
        opts([connectTimeout: '10s', requestTimeout: '45s']).httpRequestTimeout() == java.time.Duration.ofSeconds(55)

        cleanup:
        SysEnv.pop()
    }

    def 'should have no http request timeout when the request timeout is unbounded' () {
        given:
        SysEnv.push([:])

        expect:
        opts([requestTimeout: '0s']).httpRequestTimeout() == null

        cleanup:
        SysEnv.pop()
    }
}
