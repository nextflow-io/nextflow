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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerFactoryTest extends Specification {

    def 'should create an tower observer' () {

        given:
        def session = Mock(Session)
        def factory = new TowerFactory()
        
        when:
        def result = factory.create(session)[0] as TowerClient
        then:
        session.getConfig() >> [tower: [enabled: true]]
        then:
        result.endpoint == TowerClient.DEF_ENDPOINT_URL

        when:
        result = factory.create(session)[0] as TowerClient
        then:
        session.getConfig() >> [tower: [enabled: true, endpoint:'http://foo.com/api']]
        then:
        result.endpoint == 'http://foo.com/api'

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

}
