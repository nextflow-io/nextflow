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
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskPollingMonitor
import nextflow.util.Duration
import org.gridgain.grid.Grid
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.GridGain
import org.gridgain.grid.GridNode
import org.gridgain.grid.GridProjection
import org.gridgain.grid.cache.GridCache
import org.gridgain.grid.events.GridEventType
import org.gridgain.grid.spi.discovery.GridDiscoverySpiListener
import org.jetbrains.annotations.Nullable
/**
 * Creates an instance of the GridGain node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GgConnector implements GridDiscoverySpiListener {

    @Memoized
    static GgConnector create(TaskPollingMonitor monitor) {
        new GgConnector(monitor)
    }

    final private TaskPollingMonitor monitor

    final private Session session

    private GridConfiguration cfg

    private Grid grid

    private GridCache<UUID, GgClient> allSessions

    private Thread watcher

    private GgConnector(TaskPollingMonitor monitor) {
        log.debug "Create GridGain master node"
        this.monitor = monitor
        this.session = monitor.session

        initialize()
    }

    private initialize() {
        // get the configuration
        Map executorConfig = session.config.executor ?: [:]
        cfg = new GgConfigFactory('master', executorConfig).create()

        /*
         * launch the node
         */
        grid = GridGain.start(cfg);

        watcher = topologyWatcher()

        /*
         * setup the session cache
         */
        allSessions = grid.cache(GgConfigFactory.SESSIONS_CACHE)
        allSessions.put( session.uniqueId, new GgClient(session) )

        /*
         * shutdown the instance on session termination
         */
        monitor.session.onShutdown {
            if( watcher ) watcher.interrupt()
            allSessions.remove(session.uniqueId)
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

    def GridProjection getCluster() {
        def result = grid.forAttribute("ROLE", "worker")
        log.trace "Grid nodes: ${grid.nodes().size()} -- as worker: ${result.nodes().size()}"
        return result
    }

    @Override
    @Deprecated // this could be removed
    void onDiscovery(int type, long topVer, GridNode node, Collection<GridNode> topSnapshot, @Nullable Map<Long, Collection<GridNode>> topHist) {
        log.debug "Discovery event >> ${nodeEvents.get(type)?:type} [$node]"
    }

    static nodeEvents = [
            (GridEventType.EVT_NODE_FAILED): 'NODE_FAILED',
            (GridEventType.EVT_NODE_JOINED): 'NODE_JOINED',
            (GridEventType.EVT_NODE_LEFT):   'NODE_LEFT',
            (GridEventType.EVT_NODE_METRICS_UPDATED): 'EVT_NODE_METRICS_UPDATED',
            (GridEventType.EVT_NODE_RECONNECTED):     'EVT_NODE_RECONNECTED',
            (GridEventType.EVT_NODE_SEGMENTED):       'EVT_NODE_SEGMENTED'
    ]

}
