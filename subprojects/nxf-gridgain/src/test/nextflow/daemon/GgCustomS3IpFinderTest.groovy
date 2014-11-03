package nextflow.daemon

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GgCustomS3IpFinderTest extends Specification {

    def 'test custom ip finder' () {

        def result
        def finder = [:] as GgCustomS3IpFinder

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
