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
 *
 */

package io.seqera.tower.plugin

import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerFactoryTest extends Specification {

    def 'should create an tower observer' () {

        given:
        def session = Mock(Session) { getConfig() >> [:] }
        def factory = new TowerFactory(env: [TOWER_ACCESS_TOKEN: '123'])

        when:
        def client = factory.create(session)[0] as TowerClient
        then:
        session.getConfig() >> [tower: [enabled: true]]
        then:
        client.endpoint == TowerClient.DEF_ENDPOINT_URL

        when:
        client = factory.create(session)[0] as TowerClient
        then:
        session.getConfig() >> [tower: [enabled: true, endpoint:'http://foo.com/api', accessToken: 'xyz']]
        then:
        client.endpoint == 'http://foo.com/api'

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

    def 'should create with with workspace id'() {
        given:
        def session = Mock(Session)

        //
        // the workspace id is taken from the env
        //
        when:
        def factory = new TowerFactory(env: [TOWER_WORKSPACE_ID: '100'])
        def client = (TowerClient) factory.create(session)[0]
        then:
        session.getConfig() >> [tower: [enabled: true, accessToken: 'xyz']]
        and:
        client.getWorkspaceId() == '100'

        //
        // the workspace id is taken from the config
        //
        when:
        factory = new TowerFactory(env: [:])
        client = (TowerClient) factory.create(session)[0]
        then:
        session.getConfig() >> [tower: [enabled: true, workspaceId: '200', accessToken: 'xyz']]
        and:
        client.getWorkspaceId() == '200'

        //
        // the workspace id is set both in the config and the env
        // the config has the priority
        //
        when:
        factory = new TowerFactory(env: [TOWER_WORKSPACE_ID: '100'])
        client = (TowerClient) factory.create(session)[0]
        then:
        session.getConfig() >> [tower: [enabled: true, workspaceId: '200', accessToken: 'xyz']]
        and:
        client.getWorkspaceId() == '200'

        //
        // when TOWER_WORKFLOW_ID is set is a tower launch
        // then the workspace id is only taken from the env
        //
        when:
        factory = new TowerFactory(env: [TOWER_WORKSPACE_ID: '100', TOWER_WORKFLOW_ID: '111222333', TOWER_ACCESS_TOKEN: 'xyz'])
        client = (TowerClient) factory.create(session)[0]
        then:
        session.getConfig() >> [tower: [enabled: true, workspaceId: '200', accessToken: 'xyz']]
        and:
        client.getWorkspaceId() == '100'


        //
        // when enabled is false but `TOWER_WORKFLOW_ID` is provided
        // then the client should be created
        //
        when:
        factory = new TowerFactory(env: [TOWER_WORKSPACE_ID: '100', TOWER_WORKFLOW_ID: '111222333', TOWER_ACCESS_TOKEN: 'xyz'])
        client = (TowerClient) factory.create(session)[0]
        then:
        session.getConfig() >> [tower: [enabled: false]]
        and:
        client.getWorkspaceId() == '100'
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
