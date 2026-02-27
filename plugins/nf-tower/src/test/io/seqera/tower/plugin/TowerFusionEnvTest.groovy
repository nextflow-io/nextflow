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

import static com.github.tomakehurst.wiremock.client.WireMock.*

import java.time.Instant
import java.time.temporal.ChronoUnit

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import com.github.tomakehurst.wiremock.stubbing.Scenario
import com.google.gson.GsonBuilder
import io.seqera.tower.plugin.exception.UnauthorizedException
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.script.WorkflowMetadata
import nextflow.serde.gson.InstantAdapter
import spock.lang.Shared
import spock.lang.Specification
/**
 * Test cases for the TowerFusionEnv class.
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
class TowerFusionEnvTest extends Specification {

    @Shared
    WireMockServer wireMockServer

    def setupSpec() {
        wireMockServer = new WireMockServer(0)
        wireMockServer.start()
    }

    def cleanupSpec() {
        wireMockServer.stop()
    }

    def setup() {
        wireMockServer.resetAll()
        SysEnv.push([:])  // <-- ensure the system host env does not interfere
    }

    def cleanup() {
        SysEnv.pop()      // <-- restore the system host env
    }

    static String toJson(Object obj) {
        new GsonBuilder()
            .registerTypeAdapter(Instant, new InstantAdapter())
            .create()
            .toJson(obj)
    }

    def 'should return the endpoint from the config'() {
        given: 'a session'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: 'https://tower.nf'
                ]
            ]
        }

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the endpoint has the expected value'
        provider.endpoint == 'https://tower.nf'
    }

    def 'should return the endpoint from the environment'() {
        setup:
        SysEnv.push(['TOWER_API_ENDPOINT': 'https://tower.nf'])
        Global.session = Mock(Session) {
            config >> [:]
        }

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the endpoint has the expected value'
        provider.endpoint == 'https://tower.nf'

        cleanup:
        SysEnv.pop()
    }

    def 'should return the default endpoint'() {
        when: 'session config is empty'
        Global.session = Mock(Session) {
            config >> [
                tower: [:]
            ]
        }
        def provider = new TowerFusionToken()

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session config is null'
        Global.session = Mock(Session) {
            config >> null
        }

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session config is missing'
        Global.session = Mock(Session) {
            config >> [:]
        }

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session.config.tower.endpoint is not defined'
        Global.session = Mock(Session) {
            config >> [
                tower: [:]
            ]
        }

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session.config.tower.endpoint is null'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: null
                ]
            ]
        }

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session.config.tower.endpoint is empty'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: ''
                ]
            ]
        }
        provider = new TowerFusionToken()

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL

        when: 'session.config.tower.endpoint is defined as "-"'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: '-'
                ]
            ]
        }
        provider = new TowerFusionToken()

        then: 'the endpoint has the expected value'
        provider.endpoint == TowerClient.DEF_ENDPOINT_URL
    }

    def 'should return the access token from the config'() {
        given: 'a session'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    accessToken: 'abc123'
                ]
            ]
        }

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the access token has the expected value'
        provider.accessToken == 'abc123'
    }

    def 'should return the access token from the environment'() {
        setup:
        Global.session = Mock(Session) {
            config >> [:]
        }
        SysEnv.push(['TOWER_ACCESS_TOKEN': 'abc123'])

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the access token has the expected value'
        provider.accessToken == 'abc123'

        cleanup:
        SysEnv.pop()
    }

    def 'should prefer the access token from the config'() {
        setup:
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    accessToken: 'abc123'
                ]
            ]
        }
        SysEnv.push(['TOWER_ACCESS_TOKEN': 'xyz789'])

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the access token has the expected value'
        provider.accessToken == 'abc123'

        cleanup:
        SysEnv.pop()
    }

    def 'should prefer the access token from the config despite being null'() {
        setup:
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    accessToken: null
                ]
            ]
        }
        SysEnv.push(['TOWER_ACCESS_TOKEN': 'xyz789'])

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the access token has the expected value'
        provider.accessToken == null

        cleanup:
        SysEnv.pop()
    }

    def 'should prefer the access token from the environment if TOWER_WORKFLOW_ID is set'() {
        setup:
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    accessToken: 'abc123'
                ]
            ]
        }
        SysEnv.push(['TOWER_ACCESS_TOKEN' : 'xyz789', 'TOWER_WORKFLOW_ID': '123'])

        when: 'the provider is created'
        def provider = new TowerFusionToken()

        then: 'the access token has the expected value'
        provider.accessToken == 'xyz789'

        cleanup:
        SysEnv.pop()
    }

    def 'should get a license token with config'() {
        given:
        def config = [
            enabled    : true,
            endpoint   : wireMockServer.baseUrl(),
            accessToken: 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw==',
            workspaceId: '67890'
        ]
        def session = Mock(Session)
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
        session.getConfig() >> [ tower: config ]
        session.getUniqueId() >> UUID.randomUUID()
        session.getWorkflowMetadata() >> meta
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = session
        def provider = new TowerFusionToken()
        and: 'a mock endpoint at flow create'
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/trace/create?workspaceId=${config.workspaceId}"))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withBody('{"message": "", "workflowId": "1234"}')
                )
        )
        and:
        def client = TowerFactory.client(session, SysEnv.get())
        client.onFlowCreate(session)

        and: 'a mock endpoint returning a valid token'
        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo("Bearer ${config.accessToken}"))
                .withRequestBody(matchingJsonPath('$.product', equalTo("some-product")))
                .withRequestBody(matchingJsonPath('$.version', equalTo("some-version")))
                .withRequestBody(matchingJsonPath('$.workspaceId', equalTo("67890")))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withHeader('Content-Type', 'application/json')
                        .withBody('{"signedToken":"xyz789", "expiresAt":' + expirationDate + '}')
                )
        )

        when: 'a license token is requested'
        final token = provider.getLicenseToken(PRODUCT, VERSION)

        then: 'the token has the expected value'
        token == 'xyz789'

        and: 'the request is correct'
        wireMockServer.verify(1, WireMock.postRequestedFor(WireMock.urlEqualTo("/license/token/"))
            .withHeader('Authorization', WireMock.equalTo("Bearer ${config.accessToken}")))
    }

    def 'should get a license token with environment'() {
        given:
        def accessToken = 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw=='
        def workspaceId = '67890'
        SysEnv.push([
            TOWER_WORKFLOW_ID: '12345',
            TOWER_ACCESS_TOKEN: accessToken,
            TOWER_WORKSPACE_ID: workspaceId,
            TOWER_API_ENDPOINT: wireMockServer.baseUrl()
        ])
        def session = Mock(Session)
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
        session.getConfig() >> [:]
        session.getUniqueId() >> UUID.randomUUID()
        session.getWorkflowMetadata() >> meta
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = session
        def provider = new TowerFusionToken()
        and: 'a mock endpoint at flow create'
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/trace/create?workspaceId=${workspaceId}"))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withBody('{"message": "", "workflowId": "1234"}')
                )
        )
        and:
        def client = TowerFactory.client(session, SysEnv.get())
        client.onFlowCreate(session)

        and: 'a mock endpoint returning a valid token'
        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo("Bearer $accessToken"))
                .withRequestBody(matchingJsonPath('$.product', equalTo("some-product")))
                .withRequestBody(matchingJsonPath('$.version', equalTo("some-version")))
                .withRequestBody(matchingJsonPath('$.workspaceId', equalTo("${workspaceId}")))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withHeader('Content-Type', 'application/json')
                        .withBody('{"signedToken":"xyz789", "expiresAt":' + expirationDate + '}')
                )
        )

        when: 'a license token is requested'
        final token = provider.getLicenseToken(PRODUCT, VERSION)

        then: 'the token has the expected value'
        token == 'xyz789'

        and: 'the request is correct'
        wireMockServer.verify(1, WireMock.postRequestedFor(WireMock.urlEqualTo("/license/token/"))
            .withHeader('Authorization', WireMock.equalTo("Bearer ${accessToken}")))

        cleanup:
        SysEnv.pop()
    }

    def 'should refresh the auth token on 401 and retry the request'() {
        given:
        def accessToken = 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw=='
        def workspaceId = '67890'
        SysEnv.push([
            TOWER_WORKFLOW_ID: '12345',
            TOWER_ACCESS_TOKEN: accessToken,
            TOWER_REFRESH_TOKEN: 'xyz-refresh',
            TOWER_WORKSPACE_ID: workspaceId,
            TOWER_API_ENDPOINT: wireMockServer.baseUrl()
        ])
        def session = Mock(Session)
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
        session.getConfig() >> [:]
        session.getUniqueId() >> UUID.randomUUID()
        session.getWorkflowMetadata() >> meta
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = session
        def provider = new TowerFusionToken()
        and: 'a mock endpoint at flow create'
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/trace/create?workspaceId=${workspaceId}"))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withBody('{"message": "", "workflowId": "1234"}')
                )
        )
        and:
        def client = TowerFactory.client(session, SysEnv.get())
        client.onFlowCreate(session)

        and: 'prepare stubs'

        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))

        // 1️⃣ First attempt: /license/token/ fails with 401
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo("Bearer $accessToken"))
                .inScenario("Refresh flow")
                .whenScenarioStateIs(Scenario.STARTED)
                .willReturn(WireMock.aResponse().withStatus(401))
                .willSetStateTo("Token Refreshed")
        )

        // 2️⃣ Refresh token call
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/oauth/access_token"))
                .withHeader('Content-Type', equalTo('application/x-www-form-urlencoded'))
                .withRequestBody(containing('grant_type=refresh_token'))
                .withRequestBody(containing('refresh_token=xyz-refresh'))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withHeader('Set-Cookie', 'JWT=new-abc-token; Path=/; HttpOnly')
                        .withHeader('Set-Cookie', 'JWT_REFRESH_TOKEN=new-refresh-456; Path=/; HttpOnly')
                        .withBody('{"token_type":"Bearer"}')
                )
        )

        // 3️⃣ Retry: /license/token/ succeeds
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo('Bearer new-abc-token'))
                .inScenario("Refresh flow")
                .whenScenarioStateIs("Token Refreshed")
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withHeader('Content-Type', 'application/json')
                        .withBody('{"signedToken":"xyz789", "expiresAt":' + expirationDate + '}')
                )
        )

        when:
        final token = provider.getLicenseToken(PRODUCT, VERSION)

        then:
        token == 'xyz789'

        and: 'verify that refresh endpoint was called'
        wireMockServer.verify(1, WireMock.postRequestedFor(WireMock.urlEqualTo("/oauth/access_token")))

        and: 'verify both requests to license endpoint'
        wireMockServer.verify(2, WireMock.postRequestedFor(urlEqualTo("/license/token/")))

        cleanup:
        SysEnv.pop()
    }

    def 'should throw UnauthorizedException if getting a token fails with 401'() {
        given: 'a TowerFusionEnv provider'
        def config = [
            enabled    : true,
            endpoint   : wireMockServer.baseUrl(),
            accessToken: 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw==',
            workspaceId: '67890'
        ]
        def session = Mock(Session)
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
        session.getConfig() >> [ tower: config ]
        session.getUniqueId() >> UUID.randomUUID()
        session.getWorkflowMetadata() >> meta
        and:
        Global.session = session
        def provider = new TowerFusionToken()
        and: 'a mock endpoint at flow create'
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/trace/create?workspaceId=${config.workspaceId}"))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withBody('{"message": "", "workflowId": "1234"}')
                )
        )
        and:
        def client = TowerFactory.client(session, SysEnv.get())
        client.onFlowCreate(session)
        and: 'a mock endpoint returning an error'
        wireMockServer.stubFor(
            WireMock.post(WireMock.urlEqualTo("/license/token/"))
                .withHeader('Authorization', WireMock.equalTo("Bearer ${config.accessToken}"))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(401)
                        .withHeader('Content-Type', 'application/json')
                        .withBody('{"error":"Unauthorized"}')
                )
        )

        when: 'a license token is requested'
        provider.getLicenseToken('some-product', 'some-version')

        then: 'an exception is thrown'
        thrown(UnauthorizedException)
    }

    def 'should deserialize response' () {
        given:
        def ts = Instant.ofEpochSecond(1738788914)
        def json = '{"signedToken":"foo","expiresAt":"2025-02-05T20:55:14Z"}'

        when:
        def resp = TowerFusionToken.parseLicenseTokenResponse(json)
        then:
        resp.signedToken == 'foo'
        resp.expiresAt == ts
    }

}
