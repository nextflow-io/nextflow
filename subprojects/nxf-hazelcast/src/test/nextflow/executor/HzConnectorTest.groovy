/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.executor

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzConnectorTest extends  Specification {

    def testGetJoinAddress() {

        when:
        def connector = [:] as HzConnector
        connector.session = new Session( executor: [ join:'172.1.6.2' ] )
        then:
        connector.getHzJoinAddress() == ['172.1.6.2']

        when:
        connector = [:] as HzConnector
        connector.session = new Session( executor: [ join:'172.1.6.2, 172.1.6.4' ] )
        then:
        connector.getHzJoinAddress() == ['172.1.6.2','172.1.6.4']

        when:
        connector = [:] as HzConnector
        connector.session = new Session( executor: [ join:' 172.1.6.2 172.1.6.4' ] )
        then:
        connector.getHzJoinAddress() == ['172.1.6.2','172.1.6.4']

        when:
        connector = [:] as HzConnector
        connector.session = new Session( executor: [ join:' 172.1.6.2\n172.1.6.3\n\r172.1.6.4 ' ] )
        then:
        connector.getHzJoinAddress() == ['172.1.6.2','172.1.6.3','172.1.6.4']

    }
}
