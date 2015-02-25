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
import nextflow.executor.GgGridFactory
import nextflow.executor.ServiceName
import nextflow.file.ggfs.GgFileSystemProvider
import nextflow.file.ggfs.GgPath
import nextflow.file.FileHelper
import nextflow.util.KryoHelper
import nextflow.util.PathSerializer
import org.gridgain.grid.Grid
import org.weakref.s3fs.S3Path

/**
 * Launch the GridGain daemon
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ServiceName('gridgain')
class GgDaemon implements DaemonLauncher {

    Grid grid

    @Override
    void launch(Map config) {
        log.info "Configuring GridGain cluster daemon"

        /*
         * register path serializer
         */
        KryoHelper.register(GgPath, PathSerializer)
        KryoHelper.register(S3Path, PathSerializer)

        /*
         * Launch grid instance
         */
        grid = new GgGridFactory('worker', config).start()

        /*
         * configure the file system
         */
        log.debug "Configuring GridGain file system"
        FileHelper.getOrCreateFileSystemFor(GgFileSystemProvider.SCHEME, [grid: grid])
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
