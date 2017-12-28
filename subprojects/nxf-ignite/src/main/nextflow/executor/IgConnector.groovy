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

package nextflow.executor
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Session
import nextflow.cloud.CloudConfig
import nextflow.daemon.IgGridFactory
import nextflow.file.FileHelper
import nextflow.file.igfs.IgFileSystemProvider
import nextflow.file.igfs.IgPath
import nextflow.processor.TaskPollingMonitor
import nextflow.scheduler.Autoscaler
import nextflow.scheduler.Scheduler
import nextflow.scheduler.SchedulerAgent
import nextflow.util.ClusterConfig
import nextflow.util.KryoHelper
import nextflow.util.PathSerializer
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.cluster.ClusterGroup
/**
 * Creates an instance of the Ignite node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgConnector {

    @Memoized
    static IgConnector create(TaskPollingMonitor monitor) {
        new IgConnector(monitor)
    }

    final private TaskPollingMonitor monitor

    final private Session session

    private Ignite grid

    private IgniteCache<UUID, RemoteSession> allSessions

    private volatile SchedulerAgent agent

    @Delegate
    private Scheduler scheduler

    private IgConnector(TaskPollingMonitor monitor) {
        log.debug "Create Ignite master node"
        this.monitor = monitor
        this.session = monitor.session

        // initialise the connector
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

        /*
         * setup the session cache
         */
        allSessions = grid.cache(IgGridFactory.SESSIONS_CACHE)
        allSessions.put( session.uniqueId, new RemoteSession(session) )

        /*
         * shutdown the instance on session termination
         */
        final clusterConfig = new ClusterConfig(session.config.cluster as Map, Const.ROLE_MASTER, System.getenv())
        boolean shutdownCluster = clusterConfig.getAttribute('shutdownOnComplete', false) as boolean
        log.debug "Cluster shutdownOnComplete: $shutdownCluster"
        monitor.session.onShutdown {
            allSessions.remove(session.uniqueId)
            shutdown(shutdownCluster)
            // close the current instance
            grid.close()
        }


        /*
         * setup internal scheduler
         */
        def masterId = grid.cluster().localNode().id()
        scheduler = new Scheduler().init(grid, monitor)
        agent = new SchedulerAgent(grid, clusterConfig, masterId).run()

        // -- create cloud config and register the autoscaler
        final cloudEnabled = clusterConfig.isCloudCluster()
        if( cloudEnabled ) {
            final cloudConfig = CloudConfig.create(session.config)
            def autoscaler = new Autoscaler(grid, cloudConfig)
            registerAutoscaler(autoscaler)
        }
    }


    ClusterGroup getCluster() {
        def result = grid.cluster().forNodes( grid.cluster().nodes() )
        return result
    }

    /**
     * Shutdown all grid nodes
     */
    void shutdown(boolean killRemoteAgents=false) {
        log.debug "Shutting down grid nodes"
        try {
            if( killRemoteAgents ) {
                shutdownRemoteAgents()
            }
            agent.close()
            shutdownScheduler()
        }
        catch( Exception e ) {
            log.debug e.message ?: e.toString()
        }
    }

}
