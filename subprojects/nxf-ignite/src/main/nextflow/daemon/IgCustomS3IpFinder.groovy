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
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.apache.ignite.spi.IgniteSpiException
import org.apache.ignite.spi.discovery.tcp.ipfinder.s3.TcpDiscoveryS3IpFinder
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@InheritConstructors
class IgCustomS3IpFinder extends TcpDiscoveryS3IpFinder {

    private String hostIp

    {
        hostIp = System.getenv('NXF_HOST_IP')
    }

    @Override
    public void registerAddresses(Collection<InetSocketAddress> address) throws IgniteSpiException {

        Collection<InetSocketAddress> newAddress = getInetAddress(hostIp, address)
        if( newAddress != address )
            log.debug "registerAddresses $newAddress (it was: $address)"

        super.registerAddresses(newAddress)
    }

    @Override public void unregisterAddresses(Collection<InetSocketAddress> address) throws IgniteSpiException {

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
