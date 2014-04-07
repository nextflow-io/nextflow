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

import nextflow.util.DaemonConfig
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzDaemonTest extends Specification{

    def testConfigObject () {

        when:
        def daemon = [:] as HzDaemon
        daemon.config = new DaemonConfig('hazlecast',[:])
        then:
        daemon.getConfigObj().getGroupConfig().getName() == HzConst.DEFAULT_GROUP_NAME
        daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled

        when:
        daemon = [:] as HzDaemon
        daemon.config = new DaemonConfig('hazelcast',[join:'multicast', group:'pippo', port:12345])
        then:
        daemon.getConfigObj().getGroupConfig().getName() == 'pippo'
        daemon.getConfigObj().getNetworkConfig().getPort() == 12345
        daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled
        !daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().enabled

        when:
        daemon = [:] as HzDaemon
        daemon.config = new DaemonConfig('hazelcast', [join:'192.178.1.2,192.178.1.3:10 \n 192.178.1.4:20'])
        then:
        daemon.getConfigObj().getGroupConfig().getName() == HzConst.DEFAULT_GROUP_NAME
        !daemon.getConfigObj().getNetworkConfig().getJoin().getMulticastConfig().enabled
        daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().enabled
        daemon.getConfigObj().getNetworkConfig().getJoin().getTcpIpConfig().getMembers() == ['192.178.1.2','192.178.1.3:10','192.178.1.4:20']
    }


}
