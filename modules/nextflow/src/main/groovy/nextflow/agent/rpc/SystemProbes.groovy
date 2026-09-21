/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.agent.rpc

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * The production {@link Probes}: the only place in this file that touches the network, the
 * filesystem or the process environment. Every method answers {@code null}/{@code false} rather
 * than throwing, so a hostile environment degrades to an error row with a message instead of an
 * exception with a stack trace.
 */
@Slf4j
@CompileStatic
class SystemProbes implements Probes {

    /**
     * A connected {@code DatagramSocket} performs a routing-table lookup and nothing else -- no
     * packet leaves the host -- so this names the address the kernel would SOURCE traffic from,
     * which is the address a remote peer would see. On a cloud instance that is already the
     * private address the metadata service would report, which is why IMDS is never consulted
     * for the address itself.
     */
    /**
     * The IPv4 and IPv6 literals the route lookup is performed against. Nothing is sent to
     * either -- a connected {@code DatagramSocket} only consults the routing table -- so these
     * name well-known GLOBAL addresses purely to select the default route. IPv6 is tried second
     * so a dual-stack host keeps answering with its IPv4 address, while an IPv6-only fabric
     * (where the IPv4 lookup raises {@code ENETUNREACH}) still resolves.
     */
    private static final List<String> ROUTE_PROBE_ADDRESSES = ['1.1.1.1', '2606:4700:4700::1111']

    @Override
    String outboundAddress() {
        for( final target : ROUTE_PROBE_ADDRESSES ) {
            final address = routeLookup(target)
            if( address )
                return address
        }
        return null
    }

    private static String routeLookup(String target) {
        DatagramSocket socket = null
        try {
            socket = new DatagramSocket()
            socket.connect(InetAddress.getByName(target), 53)
            final local = socket.getLocalAddress()
            return local != null ? local.getHostAddress() : null
        }
        catch( Exception e ) {
            log.debug "Unable to determine the driver's outbound address towards ${target} - ${e.message}"
            return null
        }
        finally {
            socket?.close()
        }
    }

    /**
     * The interfaces a CONTAINER BRIDGE owns are skipped by name. Every Linux host with docker
     * installed carries {@code docker0} at 172.17.0.1 (plus a {@code br-*} per user-defined
     * network), and counting those as alternative driver addresses would make the multi-homed
     * warning fire on every such host -- including the plain single-NIC submit node the warning
     * is meant to distinguish from the multi-homed one.
     */
    private static final Pattern VIRTUAL_INTERFACE = Pattern.compile(/^(?:docker\d*|br-.*|podman\d*|cni-podman\d*|cni\d*|virbr\d*|veth.*|tun\d*|utun\d*)$/)

    @Override
    List<String> interfaceAddresses() {
        final List<String> result = []
        try {
            for( final nic : Collections.list(NetworkInterface.getNetworkInterfaces()) ) {
                if( !nic.isUp() || nic.isLoopback() )
                    continue
                if( VIRTUAL_INTERFACE.matcher(nic.getName()).matches() )
                    continue
                for( final addr : Collections.list(nic.getInetAddresses()) ) {
                    if( addr.isLoopbackAddress() || addr.isLinkLocalAddress() || addr.isAnyLocalAddress() )
                        continue
                    result << addr.getHostAddress()
                }
            }
        }
        catch( Exception e ) {
            log.debug "Unable to enumerate the driver's network interfaces - ${e.message}"
        }
        return result
    }

}
