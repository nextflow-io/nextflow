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

package nextflow.executor

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Session
import nextflow.cloud.CloudConfig
import nextflow.daemon.IgGridFactory
import nextflow.processor.TaskPollingMonitor
import nextflow.scheduler.Autoscaler
import nextflow.scheduler.Scheduler
import nextflow.scheduler.SchedulerAgent
import nextflow.util.ClusterConfig
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.cluster.ClusterGroup
import static nextflow.Const.ROLE_MASTER

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

        final factory = new IgGridFactory(ROLE_MASTER, session.config ?: [:])
        grid = factory.start()

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
            log.warn "Unexpected error shutting down Ignite scheduler -- ${e.message ?: e.toString()}"
        }
    }

}
