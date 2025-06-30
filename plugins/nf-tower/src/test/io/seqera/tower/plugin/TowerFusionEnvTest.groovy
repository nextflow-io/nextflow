package io.seqera.tower.plugin

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import io.seqera.tower.plugin.exception.UnauthorizedException
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.script.WorkflowMetadata
import nextflow.util.GsonHelper
import spock.lang.Shared
import spock.lang.Specification

import java.time.Instant
import java.time.temporal.ChronoUnit

import static com.github.tomakehurst.wiremock.client.WireMock.*

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

    def 'should get a license token with config'() {
        given:
        def config = [
            enabled    : true,
            endpoint   : wireMockServer.baseUrl(),
            accessToken: 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw==',
            workspaceId: '67890'
        ]
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'

        and:
        def session = Mock(Session)
        session.config >> [tower: config]
        session.getUniqueId() >> UUID.randomUUID()
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
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

        and: 'a mock endpoint returning a valid token'
        final now = Instant.now()
        final expirationDate = GsonHelper.toJson(now.plus(1, ChronoUnit.DAYS))
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
            TOWER_WORKFLOW_ID : '12345',
            TOWER_ACCESS_TOKEN: accessToken,
            TOWER_WORKSPACE_ID: workspaceId,
            TOWER_API_ENDPOINT: wireMockServer.baseUrl()
        ])
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'
        and:
        def session = Mock(Session)
        session.config >> [:]
        session.getUniqueId() >> UUID.randomUUID()
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
        session.getWorkflowMetadata() >> meta

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
        final expirationDate = GsonHelper.toJson(now.plus(1, ChronoUnit.DAYS))
        wireMockServer.stubFor(
            WireMock.post(urlEqualTo("/license/token/"))
                .withHeader('Authorization', equalTo("Bearer ${accessToken}"))
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
            .withHeader('Authorization', WireMock.equalTo("Bearer ${accessToken}")))

        cleanup:
        SysEnv.pop()
    }


    def 'should throw UnauthorizedException if getting a token fails with 401'() {
        given:
        def config = [
            enabled    : true,
            endpoint   : wireMockServer.baseUrl(),
            accessToken: 'eyJ0aWQiOiAxMTkxN30uNWQ5MGFmYWU2YjhhNmFmY2FlNjVkMTQ4ZDFhM2ZlNzlmMmNjN2I4Mw==',
            workspaceId: '67890'
        ]
        def PRODUCT = 'some-product'
        def VERSION = 'some-version'

        and:
        def session = Mock(Session)
        session.config >> [tower: config]
        session.getUniqueId() >> UUID.randomUUID()
        def meta = new WorkflowMetadata(
            session: session,
            projectName: 'the-project-name',
            repository: 'git://repo.com/foo')
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

    def 'should deserialize response'() {
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
