/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.daemon.IgGridFactory
import nextflow.file.FileHelper
import nextflow.file.igfs.IgFileSystemProvider
import nextflow.file.igfs.IgPath
import nextflow.processor.TaskPollingMonitor
import nextflow.util.ClusterConfig
import nextflow.util.Duration
import nextflow.util.KryoHelper
import nextflow.util.PathSerializer
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteCompute
import org.apache.ignite.cluster.ClusterGroup
import org.apache.ignite.cluster.ClusterNode
import org.apache.ignite.events.EventType
import org.apache.ignite.lang.IgniteRunnable
import org.apache.ignite.resources.IgniteInstanceResource
import org.apache.ignite.spi.discovery.DiscoverySpiCustomMessage
import org.apache.ignite.spi.discovery.DiscoverySpiListener
import org.jetbrains.annotations.Nullable
/**
 * Creates an instance of the Ignite node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class IgConnector implements DiscoverySpiListener {

    @Memoized
    static IgConnector create(TaskPollingMonitor monitor) {
        new IgConnector(monitor)
    }

    final private TaskPollingMonitor monitor

    final private Session session

    private Ignite grid

    private IgniteCache<UUID, RemoteSession> allSessions

    private Thread watcher

    private IgConnector(TaskPollingMonitor monitor) {
        log.debug "Create Ignite master node"
        this.monitor = monitor
        this.session = monitor.session

        initialize()
    }

    /**
     * Initialise the Ignite instance:
     * 1) register the Kryo serializer
     * 2) create an instance of the igfs file system
     * 3) create a {@link RemoteSession} object for the current {@link Session}
     * 4) create the shutdown hooks
     */
    private void initialize() {

        /*
         * Register the path serializer
         */
        KryoHelper.register(IgPath, PathSerializer)

        /*
         * access to the Ignite file system to force the instantiation of a Ignite instance
         * if it is not already available
         */
        def fs = FileHelper.getOrCreateFileSystemFor(URI.create('igfs:///'))
        grid = (fs.provider() as IgFileSystemProvider).getGrid()

        // thread to watch for topology changes and update the scheduler queue accordingly
        watcher = topologyWatcher()

        /*
         * setup the session cache
         */
        allSessions = grid.cache(IgGridFactory.SESSIONS_CACHE)
        allSessions.put( session.uniqueId, new RemoteSession(session) )

        /*
         * shutdown the instance on session termination
         */
        final config = new ClusterConfig('ignite', session.config.cluster as Map, System.getenv())
        boolean shutdownOnComplete = config.getAttribute('shutdownOnComplete', false)
        log.debug "Cluster shutdownOnComplete: $shutdownOnComplete"
        monitor.session.onShutdown {
            if( watcher ) watcher.interrupt()
            if( shutdownOnComplete )
                shutdown()
            else
                allSessions.remove(session.uniqueId)
            // close the current instance
            grid.close()
        }
    }

    /**
     * Launch a thread that will watch for added/removed cluster members increasing/decreasing accordingly
     * the monitor queue
     */
    private Thread topologyWatcher() {

        def result = Thread.start {

            float xFactor = 1.4
            int slots = cluster.metrics().getTotalCpus()
            log.debug "Cluster metrics >> slots: $slots"
            int value = (slots * xFactor) + 1
            monitor.capacitySet(value)

            def max = Duration.of('10sec').toMillis()
            while( true ) {
                try {
                    Thread.sleep(max)
                }
                catch( InterruptedException e ) {
                    break
                }
                catch( Throwable e ) {
                    log.debug("Error on topology watcher",e)
                }

                slots = cluster.metrics().getTotalCpus()
                log.debug "Cluster metrics >> cpus: $slots"
                value = (slots * xFactor) + 1
                monitor.capacitySet(value)
            }
        }

        result.setName('Topology watcher')
        return result
    }

    def ClusterGroup getCluster() {
        def result = grid.cluster().forNodes( grid.cluster().nodes() )
        return result
    }

    IgniteCompute compute() {
        grid.compute()
    }

    /**
     * Shutdown all grid nodes
     */
    def void shutdown() {
        log.debug "Shutting down grid nodes"
        try {
            grid.compute(grid.cluster().forRemotes()).broadcast(new PoisonPill())
        }
        catch( Exception e ) {
            log.debug e.message ?: e.toString()
        }
    }

    @Override
    void onDiscovery(int type,
                     long topVer,
                     ClusterNode node,
                     Collection<ClusterNode> topSnapshot,
                     @Nullable Map<Long, Collection<ClusterNode>> topHist,
                     @Nullable DiscoverySpiCustomMessage data)
    {
        log.debug "Discovery event >> ${nodeEvents.get(type)?:type} [$node]"
    }

    static nodeEvents = [
            (EventType.EVT_NODE_FAILED): 'NODE_FAILED',
            (EventType.EVT_NODE_JOINED): 'NODE_JOINED',
            (EventType.EVT_NODE_LEFT):   'NODE_LEFT',
            (EventType.EVT_NODE_METRICS_UPDATED): 'EVT_NODE_METRICS_UPDATED',
            (EventType.EVT_NODE_SEGMENTED):       'EVT_NODE_SEGMENTED'
    ]

    /**
     * Runnable task used to kill remote node on shut-down
     */
    static class PoisonPill implements IgniteRunnable {

        @IgniteInstanceResource
        private transient Ignite ignite

        @Override
        void run() {
            Thread.start { ignite.close()  }
        }
    }
}
