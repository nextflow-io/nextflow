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

        when: 'session.config.tower.endpoint is empty'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: ''
                ]
            ]
        }

        then: 'the endpoint has the expected value'

        when: 'session.config.tower.endpoint is defined as "-"'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint: '-'
                ]
            ]
        }

        then: 'the endpoint has the expected value'
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
            endpoint   : wireMockServer.baseUrl(),
            accessToken: 'abc123',
            workspaceId: '67890'
        ]
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = Mock(Session) {getConfig() >> [ tower: config ] }
        and:
        def provider = new TowerFusionToken()

        and: 'a mock endpoint returning a valid token'
        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo('Bearer abc123'))
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
            .withHeader('Authorization', WireMock.equalTo('Bearer abc123')))
    }

    def 'should get a license token with environment'() {
        given:
        SysEnv.push([
            TOWER_WORKFLOW_ID: '12345',
            TOWER_ACCESS_TOKEN: 'abc123',
            TOWER_WORKSPACE_ID: '67890',
            TOWER_API_ENDPOINT: wireMockServer.baseUrl()
        ])
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = Mock(Session) {getConfig() >> [:] }
        and:
        def provider = new TowerFusionToken()

        and: 'a mock endpoint returning a valid token'
        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo('Bearer abc123'))
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
            .withHeader('Authorization', WireMock.equalTo('Bearer abc123')))

        cleanup:
        SysEnv.pop()
    }

    def 'should refresh the auth token on 401 and retry the request'() {
        given:
        SysEnv.push([
            TOWER_WORKFLOW_ID: '12345',
            TOWER_ACCESS_TOKEN: 'abc-token',
            TOWER_REFRESH_TOKEN: 'xyz-refresh',
            TOWER_WORKSPACE_ID: '67890',
            TOWER_API_ENDPOINT: wireMockServer.baseUrl()
        ])
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        Global.session = Mock(Session) { getConfig() >> [:] }
        and:
        def provider = new TowerFusionToken()

        and: 'prepare stubs'

        final now = Instant.now()
        final expirationDate = toJson(now.plus(1, ChronoUnit.DAYS))

        // 1️⃣ First attempt: /license/token/ fails with 401
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .inScenario("Refresh flow")
                .whenScenarioStateIs(Scenario.STARTED)
                .willReturn(WireMock.aResponse().withStatus(401))
                .willSetStateTo("Token Refreshed")
        )

        // 2️⃣ Refresh token call
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/oauth/access_token"))
                .inScenario("Refresh flow")
                .whenScenarioStateIs("Token Refreshed")
                .withHeader('Content-Type', equalTo('application/x-www-form-urlencoded'))
                .withRequestBody(containing('grant_type=refresh_token'))
                .withRequestBody(containing(URLEncoder.encode('xyz-refresh', 'UTF-8')))
                .willReturn(
                    WireMock.aResponse()
                        .withStatus(200)
                        .withHeader('Set-Cookie', 'JWT=new-abc-token')
                        .withHeader('Set-Cookie', 'JWT_REFRESH_TOKEN=new-refresh-456')
                )
                .willSetStateTo("Retry Ready")
        )

        // 3️⃣ Retry: /license/token/ succeeds
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .inScenario("Refresh flow")
                .whenScenarioStateIs("Retry Ready")
                .withHeader('Authorization', equalTo('Bearer new-abc-token'))
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

        and: 'the initial request was sent with the old token'
        wireMockServer.verify(1, WireMock.postRequestedFor(urlEqualTo("/license/token/"))
            .withHeader('Authorization', equalTo('Bearer abc-token'))
        )

        and: 'the refresh request was sent with correct form data'
        wireMockServer.verify(1, WireMock.postRequestedFor(WireMock.urlEqualTo("/oauth/access_token"))
            .withHeader('Content-Type', equalTo('application/x-www-form-urlencoded'))
            .withRequestBody(containing('grant_type=refresh_token'))
            .withRequestBody(containing(URLEncoder.encode('xyz-refresh', 'UTF-8')))
        )

        and: 'the retried request was sent with the new token'
        wireMockServer.verify(1, WireMock.postRequestedFor(WireMock.urlEqualTo("/license/token/"))
            .withHeader('Authorization', equalTo('Bearer new-abc-token'))
        )

        cleanup:
        SysEnv.pop()
    }
    
    def 'should throw UnauthorizedException if getting a token fails with 401'() {
        given: 'a TowerFusionEnv provider'
        Global.session = Mock(Session) {
            config >> [
                tower: [
                    endpoint   : wireMockServer.baseUrl(),
                    accessToken: 'abc123'
                ]
            ]
        }
        and:
        def provider = new TowerFusionToken()

        and: 'a mock endpoint returning an error'
        wireMockServer.stubFor(
            WireMock.post(WireMock.urlEqualTo("/license/token/"))
                .withHeader('Authorization', WireMock.equalTo('Bearer abc123'))
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
