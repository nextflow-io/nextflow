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
