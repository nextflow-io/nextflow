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
