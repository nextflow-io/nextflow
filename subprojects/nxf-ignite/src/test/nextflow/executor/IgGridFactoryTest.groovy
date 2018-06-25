/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import nextflow.daemon.IgGridFactory
import org.apache.ignite.spi.discovery.tcp.TcpDiscoverySpi
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgGridFactoryTest extends Specification {

    def test() {
        when:
        def cfg = new IgGridFactory('master', [cluster: [ tcp:[
                localAddress:'127.1.2.3',
                localPort:8888,
                ackTimeout: '11s',
                socketTimeout: '30s',
                maxAckTimeout: 55_000,
                reconnectCount: 20,
                networkTimeout: '10s',
                joinTimeout: 100

        ]]]).config()
        def tcp = cfg.getDiscoverySpi() as TcpDiscoverySpi
        then:
        tcp.getLocalAddress() == '127.1.2.3'
        //tcp.getLocalPort() == 8888
        tcp.getAckTimeout() == 11_000       // 5 sec
        tcp.getReconnectCount() == 20       // def: 10
        tcp.getNetworkTimeout() == 10_000   // def: 5 sec
        tcp.getSocketTimeout() == 30_000   // def: 2 sec
        tcp.getMaxAckTimeout() == 55_000    // def: 600 sec
        tcp.getJoinTimeout() == 100         // def: 0

    }


}
