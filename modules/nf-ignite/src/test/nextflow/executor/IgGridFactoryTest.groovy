/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import nextflow.daemon.IgGridFactory
import org.apache.ignite.spi.discovery.tcp.TcpDiscoverySpi
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgGridFactoryTest extends Specification {

    def 'should set tcp parameters'() {
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

    def 'should set failure detection parameters'() {

        when:
        def cfg = new IgGridFactory('master', [cluster: [
                failureDetectionTimeout: '20 sec',
                clientFailureDetectionTimeout: '40 sec'
        ]]).config()
        
        then:
        cfg.getFailureDetectionTimeout() == 20_000
        cfg.getClientFailureDetectionTimeout() == 40_000

    }


}
