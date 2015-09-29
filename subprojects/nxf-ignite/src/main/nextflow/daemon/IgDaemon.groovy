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

package nextflow.daemon

import groovy.util.logging.Slf4j
import nextflow.executor.IgGridFactory
import nextflow.executor.ServiceName
import nextflow.file.FileHelper
import nextflow.file.igfs.IgFileSystemProvider
import nextflow.file.igfs.IgPath
import nextflow.util.KryoHelper
import nextflow.util.PathSerializer
import org.apache.ignite.Ignite

/**
 * Launch the GridGain daemon
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ServiceName('ignite')
class IgDaemon implements DaemonLauncher {

    Ignite grid

    @Override
    void launch(Map config) {
        log.info "Configuring Apache Ignite cluster daemon"

        /*
         * register path serializer
         */
        KryoHelper.register(IgPath, PathSerializer)

        /*
         * Launch grid instance
         */
        grid = new IgGridFactory('worker', config).start()

        /*
         * configure the file system
         */
        log.debug "Configuring Apache Ignite file system"
        FileHelper.getOrCreateFileSystemFor(IgFileSystemProvider.SCHEME, [grid: grid])
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
