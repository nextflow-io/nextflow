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
import com.amazonaws.auth.BasicAWSCredentials
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Session
import nextflow.util.DaemonConfig
import nextflow.util.Duration
import nextflow.util.FileHelper
import org.apache.commons.lang.StringUtils
import org.apache.commons.lang.reflect.MethodUtils
import org.gridgain.grid.Grid
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.GridGain
import org.gridgain.grid.cache.GridCacheAtomicityMode
import org.gridgain.grid.cache.GridCacheConfiguration
import org.gridgain.grid.cache.GridCacheDistributionMode
import org.gridgain.grid.cache.GridCacheMemoryMode
import org.gridgain.grid.cache.GridCacheMode
import org.gridgain.grid.cache.GridCacheWriteSynchronizationMode
import org.gridgain.grid.ggfs.GridGgfsConfiguration
import org.gridgain.grid.ggfs.GridGgfsGroupDataBlocksKeyMapper
import org.gridgain.grid.ggfs.GridGgfsMode
import org.gridgain.grid.logger.slf4j.GridSlf4jLogger
import org.gridgain.grid.spi.collision.jobstealing.GridJobStealingCollisionSpi
import org.gridgain.grid.spi.discovery.tcp.GridTcpDiscoverySpi
import org.gridgain.grid.spi.discovery.tcp.ipfinder.multicast.GridTcpDiscoveryMulticastIpFinder
import org.gridgain.grid.spi.discovery.tcp.ipfinder.s3.GridTcpDiscoveryS3IpFinder
import org.gridgain.grid.spi.discovery.tcp.ipfinder.sharedfs.GridTcpDiscoverySharedFsIpFinder
import org.gridgain.grid.spi.discovery.tcp.ipfinder.vm.GridTcpDiscoveryVmIpFinder
import org.gridgain.grid.spi.loadbalancing.adaptive.GridAdaptiveLoadBalancingSpi
/**
 * Grid factory class. It can be used to create a {@link GridConfiguration} or the {@link Grid} instance directly
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GgGridFactory {

    static final SESSIONS_CACHE = 'allSessions'

    static final GRID_NAME = Const.APP_NAME

    private final String role

    private final DaemonConfig config

    /**
     * Create a grid factory object for the given role and configuration params
     *
     * @param role The role for the cluster member to be configured, either {@code master} or {@code worker}
     * @param config a {@code Map} holding the configuration properties to be used
     */
    GgGridFactory( String role, Map config ) {
        assert role in ['master','worker'], "Parameter 'role' can be either 'master' or 'worker'"
        log.debug "Configuration properties for role: '$role' -- ${config}"
        this.role = role
        this.config = new DaemonConfig('gridgain', config, role == 'worker' ? System.getenv() : null )
    }

    /**
     * Creates a grid factory object by using the given {@link nextflow.Session} object.
     * @param session
     */
    GgGridFactory( Session session ) {
        this('master', session.config.executor ?: [:])
    }

    /**
     * Creates teh config object and starts a GridGain instance
     *
     * @return The {@link Grid} instance
     * @throws org.gridgain.grid.GridException If grid could not be started. This exception will be thrown
     *      also if named grid has already been started.
     */
    Grid start() {
        GridGain.start( config() )
    }

    /**
     * Main factory method, creates the {@code GridConfiguration} object
     * @return
     */
    GridConfiguration config() {

        System.setProperty('GRIDGAIN_UPDATE_NOTIFIER','false')

        if( config.getAttribute('debug') as Boolean ) {
            log.debug "Debugging ENABLED"
            System.setProperty('GRIDGAIN_DEBUG_ENABLED','true')
        }



        GridConfiguration cfg

        // -- try loading a config file
        String url = config.getAttribute('config.url')
        String file = config.getAttribute('config.file')
        if( file ) {
            log.debug "GridGain config > using config file: $file"
            cfg = GridGain.loadConfiguration(file).get1()
        }

        // -- try config by reading remote URL
        else if( url ) {
            log.debug "GridGain config > using config URL: $url"
            return GridGain.loadConfiguration(new URL(url)).get1()
        }

        // fallback of default
        else {
            cfg = new GridConfiguration()
            collisionConfig(cfg)
            discoveryConfig(cfg)
            balancingConfig(cfg)
            cacheConfig(cfg)
            ggfsConfig(cfg)
        }

        final groupName = config.getAttribute( 'group', GRID_NAME ) as String
        log.debug "GridGain config > group name: $groupName"
        cfg.setGridName(groupName)
        cfg.setUserAttributes( ROLE: role )
        cfg.setGridLogger( new GridSlf4jLogger() )

//        final addresses = config.getNetworkInterfaceAddresses()
//        if( addresses ) {
//            cfg.setLocalHost( addresses.get(0) )
//        }

        // this is not really used -- just set to avoid it complaining
        cfg.setGridGainHome( System.getProperty('user.dir') )

        return cfg
    }

    /*
     * The *session* cache
     */
    protected void cacheConfig( GridConfiguration cfg ) {

        def sessionCfg = new GridCacheConfiguration()
        sessionCfg.with {
            name = SESSIONS_CACHE
            startSize = 64
            offHeapMaxMemory = 0
        }

        /*
         * set the data cache for this ggfs
         */
        // TODO improve config setting by using a map holding default values and set all of them with an iteration
        def dataCfg = new GridCacheConfiguration()
        dataCfg.with {
            name = 'ggfs-data'
            cacheMode = GridCacheMode.PARTITIONED
            atomicityMode  = GridCacheAtomicityMode.TRANSACTIONAL
            queryIndexEnabled = false
            writeSynchronizationMode = config.getAttribute('ggfs.data.writeSynchronizationMode', GridCacheWriteSynchronizationMode.PRIMARY_SYNC) as GridCacheWriteSynchronizationMode
            distributionMode = GridCacheDistributionMode.PARTITIONED_ONLY
            affinityMapper = new GridGgfsGroupDataBlocksKeyMapper(512)
            backups = config.getAttribute('ggfs.data.backups', 0) as long
            // configure Off-heap memory
            // http://atlassian.gridgain.com/wiki/display/GG60/Off-Heap+Memory
            offHeapMaxMemory = config.getAttribute('ggfs.data.offHeapMaxMemory', 0) as long
            // When storing directly off-heap it throws an exception
            // See http://stackoverflow.com/q/23399264/395921
            memoryMode = config.getAttribute('ggfs.data.memoryMode', GridCacheMemoryMode.ONHEAP_TIERED) as GridCacheMemoryMode
        }
        cfg.setCacheConfiguration(dataCfg)

        /*
         * set the meta cache for this ggfs
         */
        def metaCfg = new GridCacheConfiguration()
        metaCfg.with {
            name = 'ggfs-meta'
            cacheMode = GridCacheMode.REPLICATED
            atomicityMode  = GridCacheAtomicityMode.TRANSACTIONAL
            queryIndexEnabled = false
            writeSynchronizationMode = GridCacheWriteSynchronizationMode.FULL_SYNC
        }


        cfg.setCacheConfiguration(sessionCfg, metaCfg, dataCfg)
    }



    /*
     * ggfs configuration
     */
    protected void ggfsConfig( GridConfiguration cfg ) {

        def ggfsCfg = new GridGgfsConfiguration()
        ggfsCfg.with {
            name = 'ggfs'
            defaultMode = GridGgfsMode.PRIMARY
            metaCacheName = 'ggfs-meta'
            dataCacheName = 'ggfs-data'
            blockSize = config.getAttribute('ggfs.blockSize', 128 * 1024) as long
            perNodeBatchSize = config.getAttribute('ggfs.perNodeBatchSize', 512) as long
            perNodeParallelBatchCount = config.getAttribute('ggfs.perNodeParallelBatchCount', 16) as long
        }
        cfg.setGgfsConfiguration(ggfsCfg)

    }


    /**
     * http://atlassian.gridgain.com/wiki/display/GG60/Job+Collision+Resolution
     * http://atlassian.gridgain.com/wiki/display/GG60/Collision+Resolution+SPI
     */
    protected collisionConfig( GridConfiguration cfg ) {

        def slots = config.getAttribute('slots', Runtime.getRuntime().availableProcessors() ) as int
        def maxActivesJobs = slots * 3
        log.debug "GridGain config > setting slots: $slots -- maxActivesJobs: $maxActivesJobs"

        def strategy = new GridJobStealingCollisionSpi()
        strategy.setActiveJobsThreshold( maxActivesJobs )
        cfg.setCollisionSpi( strategy )

    }


    /**
     * the *early* load balancing strategy
     *
     * http://www.gridgain.com/javadoc/org/gridgain/grid/spi/loadbalancing/adaptive/GridAdaptiveLoadBalancingSpi.html
     * http://atlassian.gridgain.com/wiki/display/GG60/Load+Balancing+SPI
     * http://atlassian.gridgain.com/wiki/display/GG60/Load+Balancing
     */
    protected void balancingConfig( GridConfiguration cfg ) {

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

            log.debug "GridGain config > discovery AWS bucket: $bucket; access: ${accessKey.substring(0,6)}..; ${secretKey.substring(0,6)}.."
            final finder = new GridTcpDiscoveryS3IpFinder()
            finder.setAwsCredentials( new BasicAWSCredentials(accessKey, secretKey) )
            finder.setBucketName(bucket)

            discoverCfg.setIpFinder(finder)
        }
        else if( join?.startsWith('path:') ) {
            def path = FileHelper.asPath(join.substring(5).trim())
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
        else if( join?.startsWith('ip:') ) {
            def ips = StringUtils.split(join.substring(3).trim().toString(), ", \n") as List<String>
            log.debug "GridGain config > discovery IPs: ${ips.join(', ')}"
            def finder = new GridTcpDiscoveryVmIpFinder()
            finder.setAddresses(ips)
            discoverCfg.setIpFinder(finder)

        }
        else if( join ) {
            log.warn "GrigGain config > unknown discovery method: $join"
        }

        // check some optional params
        config.getAttributesNames('tcp').each {
            checkAndSet(discoverCfg,'tcp.' + it )
        }
        cfg.setDiscoverySpi( discoverCfg )

    }

    protected void checkAndSet( def discoverCfg, String name, defValue = null ) {
        def value = config.getAttribute(name, defValue)
        if( value != null ) {
            def p = name.split(/\./)[-1]
            def x = value instanceof Duration ? value.toMillis() : value
            def n = 'set' + StringUtils.capitalize(p)
            log.debug "GridGain config > $name [$n]: $x [${x.class.name}]"
            MethodUtils.invokeMethod(discoverCfg, n, x)
        }
    }



}
