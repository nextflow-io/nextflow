/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import com.amazonaws.auth.BasicAWSCredentials
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.scheduler.JobBalancerSpi
import nextflow.scheduler.JobFailoverSpi
import nextflow.scheduler.JobSchedulerSpi
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import org.apache.commons.lang.StringUtils
import org.apache.commons.lang.reflect.MethodUtils
import org.apache.ignite.Ignite
import org.apache.ignite.Ignition
import org.apache.ignite.cache.CacheAtomicityMode
import org.apache.ignite.cache.CacheMemoryMode
import org.apache.ignite.cache.CacheMode
import org.apache.ignite.cache.CacheWriteSynchronizationMode
import org.apache.ignite.cache.eviction.lru.LruEvictionPolicy
import org.apache.ignite.configuration.CacheConfiguration
import org.apache.ignite.configuration.FileSystemConfiguration
import org.apache.ignite.configuration.IgniteConfiguration
import org.apache.ignite.igfs.IgfsGroupDataBlocksKeyMapper
import org.apache.ignite.igfs.IgfsMode
import org.apache.ignite.logger.slf4j.Slf4jLogger
import org.apache.ignite.spi.discovery.tcp.TcpDiscoverySpi
import org.apache.ignite.spi.discovery.tcp.ipfinder.multicast.TcpDiscoveryMulticastIpFinder
import org.apache.ignite.spi.discovery.tcp.ipfinder.s3.TcpDiscoveryS3IpFinder
import org.apache.ignite.spi.discovery.tcp.ipfinder.sharedfs.TcpDiscoverySharedFsIpFinder
import org.apache.ignite.spi.discovery.tcp.ipfinder.vm.TcpDiscoveryVmIpFinder
/**
 * Grid factory class. It can be used to create a {@link IgniteConfiguration} or the {@link Ignite} instance directly
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgGridFactory {

    static final public String SESSIONS_CACHE = 'allSessions'

    static final public String RESOURCE_CACHE = 'resourceCache'

    static final public String GRID_NAME = Const.APP_NAME

    static final public String NODE_ROLE = 'ROLE'

    static final public String ROLE_MASTER = 'master'

    static final public String ROLE_WORKER = 'worker'

    final private String role

    // cluster related config
    private final ClusterConfig clusterConfig

    // application config
    private final Map config

    /**
     * Create a grid factory object for the given role and configuration params
     *
     * @param role The role for the cluster member to be configured, either {@code master} or {@code worker}
     * @param config a {@code Map} holding the configuration properties to be used
     */
    IgGridFactory( String role, Map config ) {
        assert role in [ROLE_MASTER, ROLE_WORKER], "Parameter 'role' can be either `$ROLE_MASTER` or `$ROLE_WORKER`"

        final configMap = (Map)config.cluster ?: [:]
        log.debug "Configuration properties for role: '$role' -- ${configMap}"

        this.role = role
        this.clusterConfig = new ClusterConfig('ignite', configMap, System.getenv())
    }


    /**
     * Creates teh config object and starts a Ignite instance
     *
     * @return The {@link Ignite} instance
     * @throws org.apache.ignite.IgniteException If grid could not be started. This exception will be thrown
     *      also if named grid has already been started.
     */
    Ignite start() {
        Ignition.start( config() )
    }

    /**
     * Main factory method, creates the {@code IgniteConfiguration} object
     * @return
     */
    IgniteConfiguration config() {

        System.setProperty('IGNITE_UPDATE_NOTIFIER','false')
        System.setProperty('IGNITE_NO_ASCII', 'true')

        IgniteConfiguration cfg = new IgniteConfiguration()
        schedulerConfig(cfg)
        discoveryConfig(cfg)
        cacheConfig(cfg)
        fileSystemConfig(cfg)

        final groupName = clusterConfig.getAttribute( 'group', GRID_NAME ) as String
        log.debug "Apache Ignite config > group name: $groupName"
        cfg.setGridName(groupName)
        cfg.setUserAttributes( (NODE_ROLE): role )
        cfg.setGridLogger( new Slf4jLogger() )

//        final addresses = config.getNetworkInterfaceAddresses()
//        if( addresses ) {
//            cfg.setLocalHost( addresses.get(0) )
//        }

        // this is not really used -- just set to avoid it complaining
        cfg.setWorkDirectory( FileHelper.getLocalTempPath().resolve('ignite').toString() )

        return cfg
    }

    /*
     * The *session* cache
     */
    protected void cacheConfig( IgniteConfiguration cfg ) {

        def sessionCfg = new CacheConfiguration()
        sessionCfg.with {
            name = SESSIONS_CACHE
            startSize = 64
            offHeapMaxMemory = 0
        }

        /*
         * set the data cache for this ggfs
         */
        // TODO improve config setting by using a map holding default values and set all of them with an iteration
        def dataCfg = new CacheConfiguration()
        dataCfg.with {
            name = 'igfs-data'
            cacheMode = CacheMode.PARTITIONED
            evictionPolicy = new LruEvictionPolicy()
            atomicityMode  = CacheAtomicityMode.TRANSACTIONAL   // note: transactional is mandatory
            //queryIndexEnabled = false
            writeSynchronizationMode = clusterConfig.getAttribute('igfs.data.writeSynchronizationMode', CacheWriteSynchronizationMode.PRIMARY_SYNC) as CacheWriteSynchronizationMode
            //distributionMode = GridCacheDistributionMode.PARTITIONED_ONLY
            affinityMapper = new IgfsGroupDataBlocksKeyMapper(512)
            backups = clusterConfig.getAttribute('igfs.data.backups', 0) as int
            // configure Off-heap memory
            offHeapMaxMemory = clusterConfig.getAttribute('igfs.data.offHeapMaxMemory', 0) as long
            // When storing directly off-heap it throws an exception
            // See http://stackoverflow.com/q/23399264/395921
            memoryMode = clusterConfig.getAttribute('igfs.data.memoryMode', CacheMemoryMode.ONHEAP_TIERED) as CacheMemoryMode
        }
        cfg.setCacheConfiguration(dataCfg)

        /*
         * set the meta cache for this igfs
         */
        def metaCfg = new CacheConfiguration()
        metaCfg.with {
            name = 'igfs-meta'
            cacheMode = CacheMode.REPLICATED
            atomicityMode  = CacheAtomicityMode.TRANSACTIONAL   // note: transactional is mandatory
            writeSynchronizationMode = CacheWriteSynchronizationMode.PRIMARY_SYNC
            //queryIndexEnabled = false
        }

        /*
         * set scheduler resources cache
         */
        def resCfg = new CacheConfiguration()
        resCfg.with {
            name = RESOURCE_CACHE
            cacheMode = CacheMode.REPLICATED
        }

        cfg.setCacheConfiguration(sessionCfg, metaCfg, dataCfg, resCfg)

    }


    /*
     * igfs configuration
     */
    protected void fileSystemConfig( IgniteConfiguration cfg ) {

        def ggfsCfg = new FileSystemConfiguration()
        ggfsCfg.with {
            name = 'igfs'
            defaultMode = IgfsMode.PRIMARY
            metaCacheName = 'igfs-meta'
            dataCacheName = 'igfs-data'
            blockSize = clusterConfig.getAttribute('igfs.blockSize', 128 * 1024) as int
            perNodeBatchSize = clusterConfig.getAttribute('igfs.perNodeBatchSize', 512) as int
            perNodeParallelBatchCount = clusterConfig.getAttribute('igfs.perNodeParallelBatchCount', 16) as int
        }
        cfg.setFileSystemConfiguration(ggfsCfg)

    }

    protected schedulerConfig( IgniteConfiguration cfg ) {

        // -- config scheduler
        final scheduler = new JobSchedulerSpi()
        scheduler.with {
            stealingEnabled = clusterConfig.getAttribute('stealingEnabled', true) as boolean
            waitJobsThreshold = clusterConfig.getAttribute('waitJobsThreshold', DFLT_WAIT_JOBS_THRESHOLD) as int
            maximumStealingAttempts = clusterConfig.getAttribute('maxStealingAttempts', DFLT_MAX_STEALING_ATTEMPTS) as int
            maximumStealingExpireTime = clusterConfig.getAttribute('maxStealingExpireTime', DFLT_MAX_STEALING_EXPIRE_TIME) as Duration
        }
        cfg.setCollisionSpi(scheduler)

        // -- config failover
        final failover= new JobFailoverSpi()
        failover.with {
            maximumFailoverAttempts = clusterConfig.getAttribute('maxFailoverAttempts', DFLT_MAX_FAILOVER_ATTEMPTS) as int
        }
        cfg.setFailoverSpi(failover)

        // -- config load balancer
        final balancer = new JobBalancerSpi()
        cfg.setLoadBalancingSpi(balancer)

    }

    private discoveryConfig( IgniteConfiguration cfg ) {

        def discoverCfg = new TcpDiscoverySpi()

        // -- try to set the local address by using the interface configuration
        final addresses = clusterConfig.getNetworkInterfaceAddresses()
        if( addresses ) {
            final addr = addresses.get(0)
            final indx = addr.indexOf(':')
            if( indx == -1 ) {
                log.debug "Apache Ignite config > interface: $addr"
                discoverCfg.setLocalAddress(addr)
            }
            else {
                def host = addr.substring(0,indx)
                def port = addr.substring(indx+1) as Integer
                log.debug "Ignite config > interface: $host:$port"
                discoverCfg.setLocalAddress(host)
                discoverCfg.setLocalPort(port)
            }
        }

        // -- try to set the join/discovery mechanism
        def join = clusterConfig.getClusterJoin()
        if( join ) {
            if( join == 'multicast' ) {
                log.debug "Ignite config > default discovery multicast"
                discoverCfg.setIpFinder( new TcpDiscoveryMulticastIpFinder())
            }
            else if( join.startsWith('multicast:') ) {
                final finder = new TcpDiscoveryMulticastIpFinder()
                def address = join.replace('multicast:','')
                def parts = address.split(':')
                if( parts.size() ) {
                    log.debug "Ignite config > discovery multicast group: ${parts[0]}"
                    finder.setMulticastGroup(parts[0])
                }
                if( parts.size()==2 ) {
                    log.debug "Ignite config > discovery multicast port: ${parts[1]}"
                    finder.setMulticastPort(parts[1] as Integer)
                }
                discoverCfg.setIpFinder(finder)
            }

            else if( join.startsWith('s3:')) {
                def credentials = Global.getAwsCredentials(System.getenv(), config)
                if( !credentials )
                    throw new AbortOperationException("Missing AWS credentials -- Please add AWS access credentials to your environment by defining the variables AWS_ACCESS_KEY and AWS_SECRET_KEY or in your nextflow config file")

                def accessKey = credentials[0]
                def secretKey = credentials[1]
                def bucket = join.substring(3).trim()
                if( bucket.startsWith('/'))
                    bucket = bucket.substring(1)

                log.debug "Ignite config > discovery AWS bucket: $bucket; access: ${accessKey.substring(0,6)}..; ${secretKey.substring(0,6)}.."
                final finder = new TcpDiscoveryS3IpFinder()
                finder.setAwsCredentials( new BasicAWSCredentials(accessKey, secretKey) )
                finder.setBucketName(bucket)

                discoverCfg.setIpFinder(finder)
            }
            else if( join.startsWith('path:') ) {
                def path = FileHelper.asPath(join.substring(5).trim())
                if( path.exists() ) {
                    log.debug "Ignite config > discovery path: $path"
                }
                else {
                    log.debug "Ignite config > CREATING discovery path: $path"
                    path.mkdirs()
                }

                def finder = new TcpDiscoverySharedFsIpFinder()
                finder.setPath(path.toString())
                discoverCfg.setIpFinder(finder)
            }
            else if( join.startsWith('ip:') ) {
                def ips = StringUtils.split(join.substring(3).trim().toString(), ", \n") as List<String>
                log.debug "Apache Ignite config > discovery IPs: ${ips.join(', ')}"
                def finder = new TcpDiscoveryVmIpFinder()
                finder.setAddresses(ips)
                discoverCfg.setIpFinder(finder)

            }
            else {
                log.warn "Ignite config > unknown discovery method: $join"
            }
        }

        // check some optional params
        clusterConfig.getAttributesNames('tcp').each {
            checkAndSet(discoverCfg,'tcp.' + it )
        }
        cfg.setDiscoverySpi( discoverCfg )

    }

    protected void checkAndSet( def discoverCfg, String name, defValue = null ) {
        def value = clusterConfig.getAttribute(name, defValue)
        if( value != null ) {
            def p = name.split(/\./)[-1]
            def x = value instanceof Duration ? value.toMillis() : value
            def n = 'set' + StringUtils.capitalize(p)
            log.debug "Ignite config > $name [$n]: $x [${x.class.name}]"
            MethodUtils.invokeMethod(discoverCfg, n, x)
        }
    }

}
