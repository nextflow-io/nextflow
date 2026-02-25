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

import groovy.json.JsonOutput
import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import nextflow.Global
import nextflow.Session
import spock.lang.AutoCleanup
import spock.lang.Shared
import spock.lang.Specification

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 * Shared WireMock/session fixture for dataset resolver specs.
 */
abstract class DatasetWireMockSpec extends Specification {

    @Shared
    @AutoCleanup('stop')
    WireMockServer wireMock = new WireMockServer(0)

    def setupSpec() {
        wireMock.start()
        WireMock.configureFor('localhost', wireMock.port())
    }

    def setup() {
        wireMock.resetAll()
    }

    def cleanup() {
        Global.session = null
    }

    protected void mockSession(Map towerOverrides = [:]) {
        def endpoint = "http://localhost:${wireMock.port()}"
        def towerConfig = [endpoint: endpoint, accessToken: 'test-token', workspaceId: '12345'] + towerOverrides

        Global.session = Mock(Session) {
            getConfig() >> [tower: towerConfig]
        }
    }

    protected void stubDatasets(List<Map> datasets, String workspaceId = null) {
        def request = get(urlPathEqualTo('/datasets'))
        if (workspaceId)
            request = request.withQueryParam('workspaceId', equalTo(workspaceId))

        wireMock.stubFor(request.willReturn(okJson(JsonOutput.toJson([datasets: datasets]))))
    }

    protected void stubDatasetVersions(String datasetId, List<Map> versions, String workspaceId = null) {
        def request = get(urlPathEqualTo("/datasets/${datasetId}/versions"))
        if (workspaceId)
            request = request.withQueryParam('workspaceId', equalTo(workspaceId))

        wireMock.stubFor(request.willReturn(okJson(JsonOutput.toJson([versions: versions]))))
    }
}
