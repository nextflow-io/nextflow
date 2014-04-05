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
import org.gridgain.grid.Grid
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.GridGain
import org.gridgain.grid.GridProjection
import org.gridgain.grid.cache.GridCache

/**
 * Creates an instance of the GridGain node
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GgConnector {

    @Memoized
    static GgConnector create(TaskPollingMonitor monitor) {
        new GgConnector(monitor)
    }

    final private GridConfiguration cfg

    final private Grid grid

    final private TaskPollingMonitor monitor

    final private Session session

    final private GridCache<UUID, GgClient> allSessions

    private GgConnector(TaskPollingMonitor monitor) {
        log.debug "Create GridGain master node"
        this.monitor = monitor
        this.session = monitor.session

        // get the configuration
        cfg = GgConfig.create('master')

        // launch the node
        grid = GridGain.start(cfg);
        allSessions = grid.cache(GgConfig.SESSIONS_CACHE)
        allSessions.put( session.uniqueId, new GgClient(session) )

        // shutdown the instance on session termination
        monitor.session.onShutdown {
            allSessions.remove(session.uniqueId)
            grid.close()
        }
    }

    def GridProjection getCluster() {
        def result = grid.forAttribute("ROLE", "worker")
        log.trace "Grid nodes: ${grid.nodes().size()} -- as worker: ${result.nodes().size()}"
        return result
    }
}
