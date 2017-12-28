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

package nextflow.daemon

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgCustomS3IpFinderTest extends Specification {

    def 'test custom ip finder' () {

        def result
        def finder = [:] as IgCustomS3IpFinder

        when:
        result = finder.getInetAddress('1.2.3.4', null)
        then:
        result[0].getAddress().getHostAddress() == '1.2.3.4'
        result[0].getPort() == 47500

        when:
        result = finder.getInetAddress('127.10.20.255:4433', null)
        then:
        result[0].getAddress().getHostAddress() == '127.10.20.255'
        result[0].getPort() == 4433

        when:
        def addr = InetAddress.getByAddress( [0x01, 0x02, 0x03, 0x04] as byte[] )
        def fallback = [ new InetSocketAddress(addr,9977) ]
        result = finder.getInetAddress('10.20.30.40', fallback)
        then:
        result[0].getAddress().getHostAddress() == '10.20.30.40'
        result[0].getPort() == 9977

        when:
        addr = InetAddress.getByAddress( [0x01, 0x02, 0x03, 0x04] as byte[] )
        fallback = [ new InetSocketAddress(addr,9977) ]
        result = finder.getInetAddress(null, fallback)
        then:
        result[0].getAddress().getHostAddress() == '1.2.3.4'
        result[0].getPort() == 9977

    }


}
