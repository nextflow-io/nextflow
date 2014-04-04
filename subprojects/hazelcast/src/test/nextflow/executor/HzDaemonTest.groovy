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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzDaemonTest extends Specification{

    def testConfigObject () {

        when:
        def daemon = [:] as HzDaemon
        daemon.config = [:]
        then:
        daemon.getConfigObj().getGroupConfig().getName() == HzConst.DEFAULT_GROUP_NAME
        daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled

        when:
        daemon = [:] as HzDaemon
        daemon.config = [join:'multicast', group:'pippo', port:12345]
        then:
        daemon.getConfigObj().getGroupConfig().getName() == 'pippo'
        daemon.getConfigObj().getNetworkConfig().getPort() == 12345
        daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled
        !daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().enabled

        when:
        daemon = [:] as HzDaemon
        daemon.config = [join:'192.178.1.2,192.178.1.3:10 \n 192.178.1.4:20']
        then:
        daemon.getConfigObj().getGroupConfig().getName() == HzConst.DEFAULT_GROUP_NAME
        !daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled
        daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().enabled
        daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().getMembers() == ['192.178.1.2','192.178.1.3:10','192.178.1.4:20']
    }

    def testGetInterfaces() {

        setup:
        def daemon
        List<String> addresses

        when:
        daemon = [:] as HzDaemon
        addresses = daemon.getNetworkInterfaces()
        then:
        addresses == []


        when:
        daemon = [:] as HzDaemon
        daemon.config = [interface: '172.5.1.*,172.5.2.*']
        addresses = daemon.getNetworkInterfaces()
        then:
        addresses == ['172.5.1.*','172.5.2.*']

//        when:
//        daemon = [:] as HzDaemon
//        daemon.config = [interface: 'lo0']
//        addresses = daemon.getHzInterfaces()
//        then:
//        addresses == ['127.0.0.1']
    }

    def testDaemonProperty() {

        when:
        def daemon = [:] as HzDaemon
        daemon.config = [x:123, y:222, '$hazelcast': [y:333] ]
        then:
        daemon.getDaemonProperty('x') == 123
        daemon.getDaemonProperty('y') == 333
        daemon.getDaemonProperty('z') == null
        daemon.getDaemonProperty('z', 'alpha') == 'alpha'
        daemon.getDaemonProperty('z', 'alpha', [NXF_DAEMON_Z:'hola']) == 'hola'
        daemon.getDaemonProperty('p.q.z', null, [NXF_DAEMON_P_Q_Z:'hello']) == 'hello'

    }
}
