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

import nextflow.Global
import nextflow.exception.AbortOperationException
import spock.lang.Unroll

import static com.github.tomakehurst.wiremock.client.WireMock.*

/**
 * @author Edmund Miller
 */
class DatasetResolverTest extends DatasetWireMockSpec {

    def 'should throw when no session'() {
        given:
        Global.session = null

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)
    }

    @Unroll
    def 'should reject invalid dataset name: #datasetName'() {
        when:
        DatasetResolver.resolve(datasetName, null)

        then:
        thrown(IllegalArgumentException)

        where:
        datasetName << ['', null]
    }

    @Unroll
    def 'should fail dataset lookup when #scenario'() {
        given:
        mockSession()
        stubLookup.call()

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        def e = thrown(AbortOperationException)
        messageParts.each { part -> assert e.message.contains(part) }

        where:
        scenario                       | stubLookup                                                                                     | messageParts
        'dataset does not exist'       | { stubDatasets([[id: '1', name: 'other-data']]) }                                              | ['not found', 'other-data']
        'workspace has no datasets'    | { stubDatasets([]) }                                                                            | ['No datasets found in workspace']
        'api returns unauthorized'     | { wireMock.stubFor(get(urlPathEqualTo('/datasets')).willReturn(unauthorized())) }              | ['Access denied']
    }

    @Unroll
    def 'should fail version lookup when #scenario'() {
        given:
        mockSession()
        stubDatasets([[id: '42', name: 'my-data']])
        stubDatasetVersions('42', versions)

        when:
        DatasetResolver.resolve('my-data', requestedVersion)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains(expectedMessage)

        where:
        scenario                               | versions                                       | requestedVersion | expectedMessage
        'selected version has no backing URL'  | [[version: 1, url: null]]                      | null             | 'no backing storage URL'
        'requested version does not exist'     | [[version: 1, url: 's3://bucket/v1.csv']]      | '99'             | "Version '99' not found"
    }

    def 'should pass workspace ID as query param'() {
        given:
        mockSession(workspaceId: '12345')
        stubDatasets([], '12345')

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)

        and:
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets'))
            .withQueryParam('workspaceId', equalTo('12345')))
    }

    def 'should send bearer token in Authorization header'() {
        given:
        mockSession()
        stubDatasets([])

        when:
        DatasetResolver.resolve('my-data', null)

        then:
        thrown(AbortOperationException)

        and:
        wireMock.verify(getRequestedFor(urlPathEqualTo('/datasets'))
            .withHeader('Authorization', equalTo('Bearer test-token')))
    }
}
