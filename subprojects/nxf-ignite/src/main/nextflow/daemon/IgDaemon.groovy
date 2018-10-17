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

package nextflow.daemon

import groovy.util.logging.Slf4j
import nextflow.cli.CmdInfo
import nextflow.scheduler.SchedulerAgent
import nextflow.util.ServiceName
import sun.misc.Signal
import sun.misc.SignalHandler
import static nextflow.Const.ROLE_WORKER
/**
 * Launch the Ignite daemon
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ServiceName('ignite')
class IgDaemon implements DaemonLauncher {

    @Override
    void launch(Map config) {
        log.info "Configuring Apache Ignite cluster daemon"
        def info = CmdInfo.status( log.isTraceEnabled() )
        log.debug( '\n'+info )


        /*
         * Launch grid instance
         */
        def factory = new IgGridFactory(ROLE_WORKER, config)
        final grid = factory.start()

        /*
         * Scheduler agent
         */
        final agent = new SchedulerAgent(grid, factory.clusterConfig).run()

        final handler = { log.info "System terminated -- Stopping daemon.."; agent.close(true) } as SignalHandler
        Signal.handle(new Signal("INT"), handler);
        Signal.handle(new Signal("TERM"), handler);
        Signal.handle(new Signal("HUP"), handler);

    }


}
