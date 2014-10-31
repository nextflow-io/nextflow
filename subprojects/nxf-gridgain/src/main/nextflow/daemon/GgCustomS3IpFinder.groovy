package nextflow.daemon

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.gridgain.grid.spi.GridSpiException
import org.gridgain.grid.spi.discovery.tcp.ipfinder.s3.GridTcpDiscoveryS3IpFinder

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class GgCustomS3IpFinder extends GridTcpDiscoveryS3IpFinder {

    private String hostIp

    {
        hostIp = System.getenv('NXF_HOST_IP')
    }

    @Override
    public void registerAddresses(Collection<InetSocketAddress> address) throws GridSpiException {

        Collection<InetSocketAddress> newAddress = getInetAddress(hostIp, address)
        if( newAddress != address )
            log.debug "registerAddresses $newAddress (it was: $address)"

        super.registerAddresses(newAddress)
    }

    @Override public void unregisterAddresses(Collection<InetSocketAddress> address) throws GridSpiException {

        Collection<InetSocketAddress> newAddress = getInetAddress(hostIp, address)
        if( newAddress != address )
            log.debug "registerAddresses $newAddress (it was: $address)"

        super.unregisterAddresses(address)
    }


    @PackageScope
    Collection<InetSocketAddress> getInetAddress( String address, Collection<InetSocketAddress> fallback ) {
        if( !address ) {
            return fallback
        }

        def port = null
        def p = address.indexOf(':')
        if( p != -1 ) {
            port = address.substring(p+1).toInteger()
            address = address.substring(0,p)
        }
        def parts = address.split(/\./)
        if( parts.size() != 4 )
            throw new IllegalArgumentException("Not a valid IP address: $address")

        if( !port ) {
            port = fallback ? fallback[0].getPort() : 47500
        }

        def bytes = parts.collect { it.toInteger() as byte }
        def addr = InetAddress.getByAddress( bytes as byte[] )

        [ new InetSocketAddress(addr,port) ]
    }

}
