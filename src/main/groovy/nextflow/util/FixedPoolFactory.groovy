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

package nextflow.util

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.scheduler.DefaultPool
import groovyx.gpars.scheduler.Pool
import groovyx.gpars.util.PoolFactory
/**
 * Thread pool factory retuning a single instance of a fixed size thread pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FixedPoolFactory implements PoolFactory {

    private final int fNumberOfThreads

    private final boolean fDaemon

    @Lazy
    private DefaultPool fPool = {
        log.debug "Creating fixed pool size (daemon: $fDaemon, numOfThreads: $fNumberOfThreads)"
        new DefaultPool(fDaemon, fNumberOfThreads)
    } ()

    FixedPoolFactory( int size ) {
        fNumberOfThreads = size
        fDaemon = true
    }

    FixedPoolFactory( int size, boolean daemon ) {
        fNumberOfThreads = size
        fDaemon = daemon
    }

    DefaultPool getPool() { fPool }

    Pool createPool() { return fPool }

    Pool createPool(boolean daemon) {
        throw new UnsupportedOperationException()
    }

    Pool createPool(int numberOfThreads) {
        throw new UnsupportedOperationException()
    }

    Pool createPool(boolean daemon, int numberOfThreads) {
        throw new UnsupportedOperationException()
    }
}
