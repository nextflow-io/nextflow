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

import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.platform.PlatformHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerFactoryTest extends Specification {

    def 'should create a tower observer' () {
        given:
        // stub the Platform default-workspace lookup so no live API call is made when
        // no workspace is configured locally
        def factory = Spy(new TowerFactory(env: [TOWER_ACCESS_TOKEN: '123']))
        factory.defaultWorkspaceId(_ as Map) >> null

        when:
        def session = Mock(Session) { getConfig() >> [tower: [enabled: true]] }
        def observer = factory.create(session)[0] as TowerObserver
        then:
        observer.@client.endpoint == TowerClient.DEF_ENDPOINT_URL

        when:
        session = Mock(Session) { getConfig() >> [tower: [enabled: true, endpoint:'http://foo.com/api', accessToken: 'xyz']] }
        observer = factory.create(session)[0] as TowerObserver
        then:
        observer.@client.endpoint == 'http://foo.com/api'
    }

    def 'should use the Platform default workspace when none is configured locally' () {
        given:
        def factory = Spy(new TowerFactory(env: [TOWER_ACCESS_TOKEN: 'xyz']))
        def config = [tower: [enabled: true, accessToken: 'xyz']]
        def session = Mock(Session) { getConfig() >> config }

        when:
        def observer = (TowerObserver) factory.create(session)[0]

        then: 'the Platform default workspace is resolved and used'
        1 * factory.defaultWorkspaceId(_ as Map) >> '300'
        observer.getWorkspaceId() == '300'

        and: 'the resolved value is published into the session config'
        config.tower.workspaceId == '300'

        // this is what keeps Wave (registry credentials), the Fusion licence and the Seqera
        // executor scoped to the same workspace the run is reported into
        and: 'a consumer resolving the workspace afterwards observes the same value'
        PlatformHelper.getWorkspaceId(config.tower as Map, [:]) == '300'
        new TowerConfig(config.tower as Map, [:]).workspaceId == '300'
    }

    def 'should not publish anything when resolving to the personal workspace' () {
        given:
        def factory = Spy(new TowerFactory(env: [TOWER_ACCESS_TOKEN: 'xyz']))
        def config = [tower: [enabled: true, accessToken: 'xyz']]
        def session = Mock(Session) { getConfig() >> config }

        when:
        def observer = (TowerObserver) factory.create(session)[0]

        then:
        1 * factory.defaultWorkspaceId(_ as Map) >> null
        observer.getWorkspaceId() == null
        and: 'no workspace key is invented'
        !(config.tower as Map).containsKey('workspaceId')
    }

    def 'should not publish over a workspace the user configured' () {
        given:
        def factory = Spy(new TowerFactory(env: [TOWER_WORKSPACE_ID: '100', TOWER_ACCESS_TOKEN: 'xyz']))
        def config = [tower: [enabled: true, workspaceId: '200', accessToken: 'xyz']]
        def session = Mock(Session) { getConfig() >> config }

        when:
        def observer = (TowerObserver) factory.create(session)[0]

        then: 'the local setting is used and the Platform default is not queried'
        0 * factory.defaultWorkspaceId(_ as Map)
        observer.getWorkspaceId() == '200'

        and: 'the config the user wrote is left exactly as it was'
        (config.tower as Map).workspaceId == '200'
    }

    @Unroll
    def 'should fail when enabled but no access token is provided' () {
        given:
        def factory = new TowerFactory(env: [:])
        def session = Mock(Session) { getConfig() >> CONFIG }

        when:
        factory.create(session)
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('access token is required')
        e.message.contains('TOWER_ACCESS_TOKEN')

        where:
        CONFIG                          | _
        [fusion: [enabled: true]]       | _
        [tower: [enabled: true]]        | _
    }

    def 'should not create a tower observer' () {
        given:
        def session = Mock(Session)
        def factory = new TowerFactory()

        when:
        def result = factory.create(session)
        then:
        session.getConfig() >> [:]
        then:
        result == []
    }

    @Unroll
    def 'should create with workspace id from #SOURCE'() {
        given:
        def factory = Spy(new TowerFactory(env: ENV))
        def session = Mock(Session) { getConfig() >> [tower: CONFIG] }

        when:
        def observer = (TowerObserver) factory.create(session)[0]

        then: 'a local setting is available, so the Platform default is never queried'
        0 * factory.defaultWorkspaceId(_ as Map)
        observer.getWorkspaceId() == EXPECTED

        where:
        SOURCE                     | ENV                                                                                    | CONFIG                                                  | EXPECTED
        'the env'                  | [TOWER_WORKSPACE_ID: '100']                                                            | [enabled: true, accessToken: 'xyz']                     | '100'
        'the config'               | [:]                                                                                    | [enabled: true, workspaceId: '200', accessToken: 'xyz'] | '200'
        // the config has priority over the env
        'the config over the env'  | [TOWER_WORKSPACE_ID: '100']                                                            | [enabled: true, workspaceId: '200', accessToken: 'xyz'] | '200'
        // a Platform-driven run is authoritative: the workspace comes from the env only
        'a Platform-driven run'    | [TOWER_WORKSPACE_ID: '100', TOWER_WORKFLOW_ID: '111222333', TOWER_ACCESS_TOKEN: 'xyz'] | [enabled: true, workspaceId: '200', accessToken: 'xyz'] | '100'
        // the observer is created even when `enabled` is false, because TOWER_WORKFLOW_ID is set
        'a disabled Platform run'  | [TOWER_WORKSPACE_ID: '100', TOWER_WORKFLOW_ID: '111222333', TOWER_ACCESS_TOKEN: 'xyz'] | [enabled: false]                                        | '100'
    }

    @Unroll
    def 'should create tower http auth provider' () {
        given:
        def factory = new TowerFactory()
        and:
        def provider = factory.provider('https://tower.nf', 'xyz123')
        and:
        def conn = Spy(HttpURLConnection) {
            getURL() >> new URL(URL_STR)
        }

        expect:
        provider.authorize(conn) == EXPECTED
        and:
        conn.getRequestProperty('Authorization') == AUTH

        where:
        URL_STR                         | EXPECTED      | AUTH
        'http://foo.com'                | false         | null
        'https://tower.nf/'             | true          | 'Bearer xyz123'
        'https://tower.nf/this/that'    | true          | 'Bearer xyz123'
        'HTTPS://TOWER.NF/THIS/THAT'    | true          | 'Bearer xyz123'
    }

}
