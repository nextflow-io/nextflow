/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
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
