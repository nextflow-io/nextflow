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

package io.seqera.tower.plugin

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import nextflow.util.Duration

import java.net.http.HttpResponse
import java.time.Instant

import io.seqera.http.HxClient
import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerClientTest extends Specification {

    protected boolean aroundNow(value) {
        def now = Instant.now().toEpochMilli()
        value > now-1_000 && value <= now
    }

    def 'should parse response' () {
        given:
        def tower = new TowerClient()

        when:
        def resp = new TowerClient.Response(200, '{"status":"OK", "workflowId":"12345", "watchUrl": "http://foo.com/watch/12345"}')
        def result = tower.parseTowerResponse(resp)
        then:
        result.workflowId == '12345'
        result.watchUrl == 'http://foo.com/watch/12345'

        when:
        resp = new TowerClient.Response(500, '{"status":"OK", "workflowId":"12345"}')
        tower.parseTowerResponse(resp)
        then:
        thrown(Exception)
    }

    def 'should validate URL' () {
        given:
        def observer = new TowerClient()

        expect:
        observer.checkUrl('http://localhost') == 'http://localhost'
        observer.checkUrl('http://google.com') == 'http://google.com'
        observer.checkUrl('https://google.com') == 'https://google.com'
        observer.checkUrl('http://google.com:8080') == 'http://google.com:8080'
        observer.checkUrl('http://google.com:8080/') == 'http://google.com:8080'
        observer.checkUrl('http://google.com:8080/foo/bar') == 'http://google.com:8080/foo/bar'
        observer.checkUrl('http://google.com:8080/foo/bar/') == 'http://google.com:8080/foo/bar'
        observer.checkUrl('http://google.com:8080/foo/bar///') == 'http://google.com:8080/foo/bar'

        when:
        observer.checkUrl('ftp://localhost')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Only http and https are supported -- The given URL was: ftp://localhost'
    }

    def 'should get watch url' () {
        given:
        def observer = new TowerClient()
        expect:
        observer.getHostUrl(STR) == EXPECTED
        where:
        STR                             | EXPECTED
        'http://foo.com'                | 'http://foo.com'
        'http://foo.com:800/'           | 'http://foo.com:800'
        'https://foo.com:800/'          | 'https://foo.com:800'
        'http://foo.com:8000/this/that' | 'http://foo.com:8000'
    }

    def 'should get access token' () {
        when:
        def config = new TowerConfig([accessToken: 'abc'], [TOWER_ACCESS_TOKEN: 'xyz'])
        def client = new TowerClient(config, [TOWER_ACCESS_TOKEN: 'xyz'])
        then:
        // the token in the config overrides the one in the env
        client.getAccessToken() == 'abc'

        when:
        config = new TowerConfig([accessToken: 'abc'], [TOWER_ACCESS_TOKEN: 'xyz', TOWER_WORKFLOW_ID: '111222333'])
        client = new TowerClient(config, [TOWER_ACCESS_TOKEN: 'xyz', TOWER_WORKFLOW_ID: '111222333'])
        then:
        // the token from the env is taken because is a tower launch aka TOWER_WORKFLOW_ID is set
        client.getAccessToken() == 'xyz'

        when:
        config = new TowerConfig([:], [TOWER_ACCESS_TOKEN: 'xyz'])
        client = new TowerClient(config, [TOWER_ACCESS_TOKEN: 'xyz'])
        then:
        client.getAccessToken() == 'xyz'

        when:
        def c = new TowerClient()
        c.getAccessToken()
        then:
        thrown(AbortOperationException)
    }

    def 'should set the auth token' () {
        given:
        def http = Mock(HxClient.Builder)
        def client = new TowerClient()
        and:
        def SIMPLE = '4ffbf1009ebabea77db3d72efefa836dfbb71271'
        def BEARER = 'eyJ0aWQiOiA1fS5jZmM1YjVhOThjZjM2MTk1NjBjZWU1YmMwODUxYzA1ZjkzMDdmN2Iz'

        when:
        client.setupClientAuth(http, SIMPLE)
        then:
        1 * http.basicAuth('@token:' + SIMPLE) >> http

        when:
        client.setupClientAuth(http, SIMPLE)
        then:
        1 * http.basicAuth('@token:' + SIMPLE) >> http

        when:
        client.setupClientAuth(http, BEARER)
        then:
        1 * http.bearerToken(BEARER) >> http
        1 * http.refreshToken(_) >> http
        1 * http.refreshTokenUrl(_) >> http
    }

    def 'should get trace endpoint' () {
        given:
        def client = new TowerClient()
        client.@endpoint = TowerClient.DEF_ENDPOINT_URL

        expect:
        client.getUrlTraceCreate(null) == 'https://api.cloud.seqera.io/trace/create'
        client.getUrlTraceBegin(null, '12345') == 'https://api.cloud.seqera.io/trace/12345/begin'
        client.getUrlTraceProgress(null, '12345') == 'https://api.cloud.seqera.io/trace/12345/progress'
        client.getUrlTraceHeartbeat(null, '12345') == 'https://api.cloud.seqera.io/trace/12345/heartbeat'
        client.getUrlTraceComplete(null, '12345') == 'https://api.cloud.seqera.io/trace/12345/complete'
    }

    def 'should get trace endpoint with workspace' () {
        given:
        def client = new TowerClient()
        client.@endpoint = TowerClient.DEF_ENDPOINT_URL

        expect:
        client.getUrlTraceCreate('300') == 'https://api.cloud.seqera.io/trace/create?workspaceId=300'
        client.getUrlTraceBegin('300', '12345') == 'https://api.cloud.seqera.io/trace/12345/begin?workspaceId=300'
        client.getUrlTraceProgress('300', '12345') == 'https://api.cloud.seqera.io/trace/12345/progress?workspaceId=300'
        client.getUrlTraceHeartbeat('300', '12345') == 'https://api.cloud.seqera.io/trace/12345/heartbeat?workspaceId=300'
        client.getUrlTraceComplete('300', '12345') == 'https://api.cloud.seqera.io/trace/12345/complete?workspaceId=300'
    }

    def 'should load schema col len' () {
        given:
        def tower = new TowerClient()

        when:
        def schema = tower.loadSchema()
        then:
        schema.get('workflow.start')  == null
        schema.get('workflow.profile') == 100
        schema.get('workflow.projectDir') == 255
    }

    def 'should handle HTTP request with content'() {
        given: 'a TowerClient'
        def tower = new TowerClient()
        def content = '{"test": "data"}'
        def request = tower.makeRequest('http://example.com/test', content, 'POST')

        expect: 'the request should be created with the content'
        request != null
        request.method() == 'POST'
        request.uri().toString() == 'http://example.com/test'
    }

    def 'should send http message' () {
        given:
        def client = Mock(HxClient)
        def tower = new TowerClient()
        tower.@httpClient = client

        when:
        def resp = tower.sendHttpMessage('http://foo.com', [foo: 'bar'], 'POST')
        then:
        1 * client.sendAsString(_) >> Mock(HttpResponse) { statusCode() >> 200; body() >> '{}' }
        and:
        !resp.error
        resp.code == 200
    }

    def 'should return error response on http request timeout' () {
        given: 'a WireMock server that hangs for 5 seconds'
        def wireMock = new WireMockServer(0)
        wireMock.start()
        wireMock.stubFor(
            WireMock.post(WireMock.anyUrl())
                .willReturn(WireMock.aResponse()
                    .withFixedDelay(5_000)
                    .withStatus(200)
                    .withBody('{}'))
        )

        and: 'a TowerClient whose requests carry a 200ms timeout'
        TowerConfig config = Mock(TowerConfig) {
            getHttpReadTimeout() >> Duration.of('200 ms')
            getHttpConnectTimeout() >> Duration.of('5 s')
            getEndpoint() >> wireMock.baseUrl()
            getAccessToken() >> 'token'
        }
        TowerClient client = new TowerClient(config, [:])

        when:
        def response = client.sendHttpMessage("${wireMock.baseUrl()}/trace/create", [runName: 'test'], 'POST')

        then: 'a timeout produces an error response with code 0'
        response.code == 0
        response.message.contains('Unable to connect')

        cleanup:
        wireMock.stop()
    }


}
