/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import java.nio.file.Paths

import com.amazonaws.auth.BasicAWSCredentials
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.util.DaemonConfig
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.cache.GridCacheConfiguration
import org.gridgain.grid.logger.slf4j.GridSlf4jLogger
import org.gridgain.grid.spi.discovery.tcp.GridTcpDiscoverySpi
import org.gridgain.grid.spi.discovery.tcp.ipfinder.multicast.GridTcpDiscoveryMulticastIpFinder
import org.gridgain.grid.spi.discovery.tcp.ipfinder.s3.GridTcpDiscoveryS3IpFinder
import org.gridgain.grid.spi.discovery.tcp.ipfinder.sharedfs.GridTcpDiscoverySharedFsIpFinder
import org.gridgain.grid.spi.loadbalancing.adaptive.GridAdaptiveLoadBalancingSpi
/**
 * Creates a the GridGain configuration object common to master and worker nodes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GgConfigFactory {

    static final SESSIONS_CACHE = 'allSessions'

    static final GRID_NAME = Const.APP_NAME

    private final String role

    private final DaemonConfig config

    /**
     * @param role The role for the cluster member to be configured, either {@code master} or {@code worker}
     * @param config a {@code Map} holding the configuration properties to be used
     */
    GgConfigFactory( String role, Map config ) {
        assert role in ['master','worker'], "Parameter 'role' can be either 'master' or 'worker'"
        log.debug "Configuration properties for role: '$role' -- ${config}"
        this.role = role
        this.config = new DaemonConfig('gridgain', config, role == 'worker' ? System.getenv() : null )
    }

    /**
     * Main factory method, creates the {@code GridConfiguration} object
     * @return
     */
    GridConfiguration create( ) {

        // configure
        System.setProperty('GRIDGAIN_DEBUG_ENABLED','true')
        System.setProperty('GRIDGAIN_UPDATE_NOTIFIER','false')

        def cfg = new GridConfiguration()
        cfg.setGridName( GRID_NAME )
        cfg.setUserAttributes( ROLE: role )
        cfg.setGridLogger( new GridSlf4jLogger() )

        discoveryConfig(cfg)
        loadBalancingConfig(cfg)
        cacheConfig(cfg)

        return cfg
    }

    /*
     * The *session* cache
     */
    private void cacheConfig( GridConfiguration cfg ) {

        def sessionCfg = new GridCacheConfiguration()
        sessionCfg.setName(SESSIONS_CACHE)
        sessionCfg.setStartSize(64)
        cfg.setCacheConfiguration(sessionCfg)
    }


    /**
     * the *early* load balancing strategy
     *
     * http://www.gridgain.com/javadoc/org/gridgain/grid/spi/loadbalancing/adaptive/GridAdaptiveLoadBalancingSpi.html
     * http://atlassian.gridgain.com/wiki/display/GG60/Load+Balancing+SPI
     * http://atlassian.gridgain.com/wiki/display/GG60/Load+Balancing
     */
    static private void loadBalancingConfig( GridConfiguration cfg ) {

        cfg.setLoadBalancingSpi( new GridAdaptiveLoadBalancingSpi() )
    }

    /*
     * http://gridgain.com/javadoc/org/gridgain/grid/spi/discovery/tcp/GridTcpDiscoverySpi.html
     * http://atlassian.gridgain.com/wiki/display/GG60/Discovery+SPI
     */
    private discoveryConfig( GridConfiguration cfg ) {

        def discoverCfg = new GridTcpDiscoverySpi()

        // -- try to set the local address by using the interface configuration
        final addresses = config.getNetworkInterfaceAddresses()
        if( addresses ) {
            final addr = addresses.get(0)
            final indx = addr.indexOf(':')
            if( indx == -1 ) {
                log.debug "GridGain config > interface: $addr"
                discoverCfg.setLocalAddress(addr)
            }
            else {
                def host = addr.substring(0,indx)
                def port = addr.substring(indx+1) as Integer
                log.debug "GridGain config > interface: $host:$port"
                discoverCfg.setLocalAddress(host)
                discoverCfg.setLocalPort(port)
            }
        }

        // -- try to set the join/discovery mechanism
        def join = config.getAttribute('join') as String
        if( join == 'multicast' ) {
            log.debug "GridGain config > default discovery multicast"
            discoverCfg.setIpFinder( new GridTcpDiscoveryMulticastIpFinder())
        }
        else if( join?.startsWith('multicast:') ) {
            final finder = new GridTcpDiscoveryMulticastIpFinder()
            def address = join.replace('multicast:','')
            def parts = address.split(':')
            if( parts.size() ) {
                log.debug "GridGain config > discovery multicast address: ${parts[0]}"
                finder.setMulticastGroup(parts[0])
            }
            if( parts.size()==2 ) {
                log.debug "GridGain config > discovery multicast port: ${parts[1]}"
                finder.setMulticastPort(parts[1] as Integer)
            }
            discoverCfg.setIpFinder(finder)
        }

        else if( join?.startsWith('s3:')) {
            String accessKey = config.getAttribute('awsAccessKey')
            String secretKey = config.getAttribute('awsSecretKey')
            if( !accessKey ) accessKey = System.getenv('AWS_ACCESS_KEY')
            if( !secretKey ) secretKey = System.getenv('AWS_SECRET_KEY')

            def bucket = join.substring(3).trim()
            if( bucket.startsWith('/'))
                bucket = bucket.substring(1)

            log.debug "GridGain config > discovery AWS bucket: $bucket; access: ${accessKey.substring(6)}..; ${secretKey.substring(6)}.."
            final finder = new GridTcpDiscoveryS3IpFinder()
            finder.setAwsCredentials( new BasicAWSCredentials(accessKey, secretKey) )
            finder.setBucketName(bucket)

            discoverCfg.setIpFinder(finder)
        }
        else if( join?.startsWith('file:') ) {
            def path = Paths.get(join.substring(5).trim())
            if( path.exists() ) {
                log.debug "GridGain config > discovery path: $path"
            }
            else {
                log.debug "GridGain config > CREATING discovery path: $path"
                path.mkdirs()
            }

            def finder = new GridTcpDiscoverySharedFsIpFinder()
            finder.setPath(path.toString())
            discoverCfg.setIpFinder(finder)
        }
        else if( join ) {
            log.warn "GrigGain config > unknown discovery method: $join"
        }

        cfg.setDiscoverySpi( discoverCfg )
    }




}
