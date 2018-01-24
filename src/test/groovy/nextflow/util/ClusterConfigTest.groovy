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

package nextflow.util

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ClusterConfigTest extends Specification {
    
    def testGetAttribute() {

        when:
        def cfg = new ClusterConfig([x:123, y:222, 'master': [y:333] ], 'master')
        then:
        cfg.getAttribute('x') == 123
        cfg.getAttribute('y') == 333
        cfg.getAttribute('z') == null
        cfg.getAttribute('z', 'alpha') == 'alpha'

        when:
        def env = [NXF_CLUSTER_Z:'hola', NXF_CLUSTER_P_Q_Z:'hello']
        cfg = new ClusterConfig([x:123, y:222, 'master': [y:333] ], 'master',  env )
        then:
        cfg.getAttribute('z', 'alpha') == 'hola'
        cfg.getAttribute('p.q.z', null) == 'hello'
        cfg.getAttribute('y') == 333

    }

    def testGetAttributeWithType() {

        given:
        def cfg = new ClusterConfig([alpha:123, beta:'222', gamma: 'false'], 'myDaemon' )

        when:
        String str = cfg.getAttribute('alpha')
        then:
        str == '123'

        when:
        Integer num = cfg.getAttribute('alpha')
        then:
        num == 123

        when:
        Integer num_beta = cfg.getAttribute('beta')
        then:
        num_beta == 222

        when:
        boolean flag = cfg.getAttribute('gamma')
        then:
        flag == false
    }

    def testGetAttributeFromEnvironmentVariable() {

        given:
        def cfg = new ClusterConfig(null, 'ignite', [NXF_CLUSTER_SHUTDOWN: 'true'] )
        expect:
        (boolean)cfg.getAttribute('shutdown') == true
        cfg.getAttribute('something') == null
    }


    def testGetInterfaces() {

        setup:
        def config
        List<String> addresses

        when:
        config = new ClusterConfig([:])
        addresses = config.getNetworkInterfaceAddresses()
        then:
        addresses == []


        when:
        config = new ClusterConfig([interface: '172.5.1.*,172.5.2.*'])
        addresses = config.getNetworkInterfaceAddresses()
        then:
        addresses == ['172.5.1.*','172.5.2.*']

    }

    def testGetNestedNames() {

        when:
        def cfg = new ClusterConfig([x:1, y:2, tcp: [alpha: 'a', beta: 'b', gamma: [uno:1, due: 2]] ])
        then:
        cfg.getAttributeNames() == ['x','y','tcp'] as Set
        cfg.getAttributeNames('tcp') == ['alpha','beta', 'gamma'] as Set
        cfg.getAttributeNames('tcp.gamma') == ['uno','due'] as Set

    }


    def testNestedValues() {

        when:
        def cfg = new ClusterConfig([x:1, y:2, tcp: [alpha: 'a', beta: 'b', gamma: [uno:1, due: 2]] ])
        then:
        cfg.getAttribute('x') == 1
        cfg.getAttribute('tcp.alpha') == 'a'
        cfg.getAttribute('tcp.beta') == 'b'
        cfg.getAttribute('tcp.gamma') == [ uno:1, due: 2]
        cfg.getAttribute('tcp.gamma.uno') == 1
        cfg.getAttribute('tcp.gamma.due') == 2

    }

    enum FooEnum { ALPHA, BETA, DELTA }

    def testEnumValues() {

        when:
        def cfg = new ClusterConfig([x:1, ggfs:[data: [memoryMode: 'ALPHA']] ])
        then:
        cfg.getAttribute('ggfs.data.memoryMode') == 'ALPHA'
        cfg.getAttribute('ggfs.data.memoryMode') as FooEnum == FooEnum.ALPHA
        cfg.getAttribute('ggfs.meta.memoryMode', FooEnum.BETA) as FooEnum == FooEnum.BETA
        cfg.getAttribute('ggfs.meta.xxx') as FooEnum == null

    }

    @Unroll
    def 'should return multicast group address for groupId: #groupId' () {
        expect:
        ClusterConfig.seedToMulticastAddress(groupId) == address
        where:
        groupId     | address
        '1'         | 'multicast:228.0.0.1'
        '255'       | 'multicast:228.0.0.255'
        '256'       | 'multicast:228.0.1.0'
        '65535'     | 'multicast:228.0.255.255'
        '65536'     | 'multicast:228.1.0.0'
        '65792'     | 'multicast:228.1.1.0'
        '65793'     | 'multicast:228.1.1.1'
        '66052'     | 'multicast:228.1.2.4'
        '16777215'  | 'multicast:228.255.255.255'
        '16777216'  | null

    }

    def 'should return cluster join' () {
        when:
        def cfg = new ClusterConfig([join: 'path:/some/dir'])
        then:
        cfg.getClusterJoin() == 'path:/some/dir'

        when:
        cfg = new ClusterConfig(null, 'ignite', [NXF_CLUSTER_JOIN: 's3://bucket'])
        then:
        cfg.getClusterJoin() == 's3://bucket'

        // fallback on NXF_CLUSTER_SEED
        when:
        cfg = new ClusterConfig(null, 'ignite', [NXF_CLUSTER_SEED: '66052'])
        then:
        cfg.getClusterJoin() == 'multicast:228.1.2.4'
    }


}
