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
import groovy.util.logging.Slf4j
import nextflow.daemon.DaemonLauncher
import org.gridgain.grid.Grid
import org.gridgain.grid.GridConfiguration
import org.gridgain.grid.GridGain

/**
 * Launch the GridGain daemon
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class GgDaemon implements DaemonLauncher {

    Grid grid

    GridConfiguration cfg

    @Override
    void launch(Map config) {
        log.info "Configuring GridGain cluster daemon"
        run()
    }

    void run() {
        cfg = GgConfig.create('worker')
        grid = GridGain.start(cfg);
    }



//    @Slf4j
//    class DaemonLifecycleBean implements GridLifecycleBean {
//
//        @Override
//        void onLifecycleEvent(GridLifecycleEventType event) throws GridException {
//            if( event == GridLifecycleEventType.BEFORE_GRID_START ) {
//                onBeforeStart()
//            }
//            else if( event == GridLifecycleEventType.AFTER_GRID_START ) {
//                onAfterStart()
//            }
//            else if( event == GridLifecycleEventType.BEFORE_GRID_STOP ) {
//                onBeforeStop()
//            }
//            else if( event == GridLifecycleEventType.AFTER_GRID_STOP ) {
//                onAfterStop()
//            }
//        }
//
//        def onBeforeStart() {
//            log.debug "Grid > before start"
//        }
//
//        def onAfterStart() {
//            log.debug "Grid > after start"
//        }
//
//        def onBeforeStop() {
//            log.debug "Grid > before stop"
//        }
//
//        def onAfterStop() {
//            log.debug "Grid > after stop"
//        }
//    }
}
