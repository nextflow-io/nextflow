package nextflow.executor

import org.gridgain.grid.spi.discovery.tcp.GridTcpDiscoverySpi
import spock.lang.Specification

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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GgConfigFactoryTest extends Specification {

    def test() {
        when:
        def cfg = new GgConfigFactory('master', [tcp:[
                localAddress:'127.1.2.3',
                localPort:8888,
                ackTimeout: '11s',
                socketTimeout: '30s',
                maxAckTimeout: 55_000,
                heartbeatFrequency: '3sec',
                maxMissedHeartbeats: 3,
                reconnectCount: 20,
                networkTimeout: '10s',
                joinTimeout: 100,

            ], ]).create()
        def tcp = cfg.getDiscoverySpi() as GridTcpDiscoverySpi
        then:
        tcp.getLocalAddress() == '127.1.2.3'
        //tcp.getLocalPort() == 8888
        tcp.getAckTimeout() == 11_000       // 5 sec
        tcp.getHeartbeatFrequency() == 3_000    // 2 sec
        tcp.getMaxMissedHeartbeats() == 3       // 1
        tcp.getReconnectCount() == 20       // def: 10
        tcp.getNetworkTimeout() == 10_000   // def: 5 sec
        tcp.getSocketTimeout() == 30_000   // def: 2 sec
        tcp.getMaxAckTimeout() == 55_000    // def: 600 sec
        tcp.getJoinTimeout() == 100         // def: 0

    }


}
