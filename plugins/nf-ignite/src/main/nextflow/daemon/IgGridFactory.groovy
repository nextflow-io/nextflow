/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import static nextflow.Const.*

import java.util.logging.Level
import java.util.logging.Logger

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.file.FileHelper
import nextflow.ignite.IgnitePlugin
import nextflow.scheduler.Protocol
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import org.apache.commons.lang.StringUtils
import org.apache.commons.lang.reflect.MethodUtils
import org.apache.ignite.Ignite
import org.apache.ignite.Ignition
import org.apache.ignite.cache.CacheMode
import org.apache.ignite.configuration.CacheConfiguration
import org.apache.ignite.configuration.IgniteConfiguration
import org.apache.ignite.logger.slf4j.Slf4jLogger
import org.apache.ignite.spi.discovery.tcp.TcpDiscoverySpi
import org.apache.ignite.spi.discovery.tcp.ipfinder.multicast.TcpDiscoveryMulticastIpFinder
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

    static final public String GRID_NAME = Const.APP_NAME

    static final public String NODE_ROLE = 'ROLE'

    final private String role

    // cluster related config
    private final ClusterConfig clusterConfig

    // application config
    private final Map config

    static private IgGridFactory singleton

    static Ignite ignite() { Ignition.ignite(GRID_NAME) }

    static IgGridFactory instance() { singleton }

    ClusterConfig getClusterConfig() { clusterConfig }

    /**
     * Create a grid factory object for the given role and configuration params
     *
     * @param role The role for the cluster member to be configured, either {@code master} or {@code worker}
     * @param config a {@code Map} holding the configuration properties to be used
     */
    IgGridFactory( String role, Map config ) {
        assert role in [ROLE_MASTER, ROLE_WORKER], "Parameter 'role' can be either `$ROLE_MASTER` or `$ROLE_WORKER`"
        singleton = this

        final configMap = (Map)config.cluster ?: [:]
        log.debug "Configuration properties for role: '$role' -- ${configMap}"

        this.role = role
        this.clusterConfig = new ClusterConfig(configMap, role, System.getenv())
    }


    /**
     * Creates teh config object and starts a Ignite instance
     *
     * @return The {@link Ignite} instance
     * @throws org.apache.ignite.IgniteException If grid could not be started. This exception will be thrown
     *      also if named grid has already been started.
     */
    Ignite start() {
        final oldLdr = Thread.currentThread().getContextClassLoader()
        try {
            final plgLdr = IgnitePlugin.class.getClassLoader()
            log.debug "+ Setting context class loader to=$plgLdr - previous=$oldLdr"
            Thread.currentThread().setContextClassLoader(plgLdr)
            Ignition.start( config() )
        }
        finally {
            Thread.currentThread().setContextClassLoader(oldLdr)
        }
    }

    /**
     * Main factory method, creates the {@code IgniteConfiguration} object
     * @return
     */
    IgniteConfiguration config() {
        // required by
        // https://issues.apache.org/jira/browse/IGNITE-8899
        // https://issues.apache.org/jira/browse/IGNITE-8426
        Logger.getLogger('').setLevel(Level.OFF)

        System.setProperty('IGNITE_UPDATE_NOTIFIER','false')
        System.setProperty('IGNITE_NO_ASCII', 'true')
        System.setProperty('IGNITE_NO_SHUTDOWN_HOOK', 'true')
        System.setProperty('IGNITE_QUIET', 'false')

        IgniteConfiguration cfg = new IgniteConfiguration()
        // set the classloader to allow loading classes from the plugin path
        log.debug "+ Setting Ignite plugin classloader=${IgnitePlugin.class.getClassLoader()}"
        cfg.setClassLoader(IgnitePlugin.class.getClassLoader())

        discoveryConfig(cfg)
        cacheConfig(cfg)

        final groupName = clusterConfig.getAttribute( 'group', GRID_NAME ) as String
        log.debug "Apache Ignite config > group name: $groupName"
        cfg.setIgniteInstanceName(groupName)
        cfg.setUserAttributes( (NODE_ROLE): role )
        cfg.setGridLogger( new Slf4jLogger() )

        def freq = clusterConfig.getAttribute('metricsLogFrequency', Duration.of('5 min')) as Duration
        cfg.setMetricsLogFrequency(freq.toMillis())

        // this is not really used -- just set to avoid it complaining
        cfg.setWorkDirectory( FileHelper.getLocalTempPath().resolve('ignite').toString() )

        def timeout = clusterConfig.getAttribute('failureDetectionTimeout') as Duration
        if( timeout )
            cfg.setFailureDetectionTimeout(timeout.millis)

        timeout = clusterConfig.getAttribute('clientFailureDetectionTimeout') as Duration
        if( timeout )
            cfg.setClientFailureDetectionTimeout(timeout.millis)

        log.debug "Apache Ignite config > $clusterConfig"
        
        return cfg
    }

    /*
     * The *session* cache
     */
    protected void cacheConfig( IgniteConfiguration cfg ) {

        List<CacheConfiguration> configs = []

        configs << new CacheConfiguration()
                .setName(SESSIONS_CACHE)

        /*
         * set scheduler resources cache
         */
        configs << new CacheConfiguration()
                .setName(Protocol.PENDING_TASKS_CACHE)
                .setCacheMode(CacheMode.REPLICATED)

        cfg.setCacheConfiguration( configs as CacheConfiguration[] )
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
        clusterConfig.getAttributeNames('tcp').each {
            checkAndSet(discoverCfg,'tcp.' + it )
        }
        cfg.setDiscoverySpi( discoverCfg )

    }

    private String getLocalAddress() {
        try {
            return InetAddress.getLocalHost().getHostAddress()
        }
        catch( IOException e ) {
            log.debug "Oops.. Cannot find local address", e
            return null
        }

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
