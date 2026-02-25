/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.tower.plugin.dataset

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import spock.lang.AutoCleanup
import spock.lang.Specification

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 * @author Edmund Miller
 */
class DatasetResolverTest extends Specification {

    @AutoCleanup('stop')
    WireMockServer wireMock

    def setup() {
        wireMock = new WireMockServer(0)
        wireMock.start()
        WireMock.configureFor(wireMock.port())
    }

    def cleanup() {
        Global.session = null
    }

    private void mockSession(Map extra = [:]) {
        def endpoint = "http://localhost:${wireMock.port()}"
        def config = [tower: [endpoint: endpoint, accessToken: 'test-token', workspaceId: '12345'] + extra]
        Global.session = Mock(Session) {
            getConfig() >> config
        }
    }

    def 'should throw when no session'() {
        given:
        Global.session = null

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)
    }

    def 'should throw when dataset name is empty'() {
        when:
        DatasetResolver.resolve('', null)

        then:
        thrown(IllegalArgumentException)
    }

    def 'should throw when dataset not found'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "1", "name": "other-data"}]}')))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("not found")
        e.message.contains("other-data")
    }

    def 'should throw when no datasets in workspace'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": []}')))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)
    }

    def 'should throw when API returns 401'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(unauthorized()))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Access denied")
    }

    def 'should throw when version has no backing URL'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "42", "name": "my-data"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/42/versions'))
            .willReturn(okJson('{"versions": [{"version": 1, "url": null}]}')))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("no backing storage URL")
    }

    def 'should throw when specific version not found'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": [{"id": "42", "name": "my-data"}]}')))
        wireMock.stubFor(get(urlPathEqualTo('/datasets/42/versions'))
            .willReturn(okJson('{"versions": [{"version": 1, "url": "s3://bucket/v1.csv"}]}')))

        when:
        DatasetResolver.resolve('my-data', '99')

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Version '99' not found")
    }

    def 'should pass workspace ID as query param'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .withQueryParam('workspaceId', equalTo('12345'))
            .willReturn(okJson('{"datasets": []}')))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        // will throw because no datasets, but the important thing
        // is the request was made with correct query param
        thrown(AbortOperationException)

        and:
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets'))
            .withQueryParam('workspaceId', equalTo('12345')))
    }

    def 'should send bearer token in Authorization header'() {
        given:
        mockSession()

        wireMock.stubFor(get(urlPathEqualTo('/datasets'))
            .willReturn(okJson('{"datasets": []}')))

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)

        and:
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets'))
            .withHeader('Authorization', equalTo('Bearer test-token')))
    }
}
