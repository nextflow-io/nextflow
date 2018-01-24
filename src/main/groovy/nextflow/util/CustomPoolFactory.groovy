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

package nextflow.util
import groovy.transform.CompileStatic
import groovyx.gpars.scheduler.Pool
import groovyx.gpars.scheduler.ResizeablePool
import groovyx.gpars.util.PoolFactory
/**
 * A configurable thread pool factory
 *
 * See http://niklasschlimm.blogspot.com.es/2012/03/threading-stories-about-robust-thread.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CustomPoolFactory implements PoolFactory {

    public final String NXF_POOL_TYPE = 'nxf.pool.type'

    public final String NXF_MAX_THREADS = 'nxf.pool.maxThreads'

    public final String NXF_QUEUE_SIZE = 'nxf.pool.queueSize'

    @Override
    Pool createPool() {

        def type = property(NXF_POOL_TYPE, 'default')
        switch (type) {
            case 'default':
                return new ResizeablePool(true, 1)

            case 'sync':
                int cpus = Runtime.runtime.availableProcessors()
                int size = property(NXF_MAX_THREADS, cpus+1) as int
                return CustomThreadPool.synchronousPool(size)

            case 'bound':
                int cpus = Runtime.runtime.availableProcessors()
                int size = property(NXF_MAX_THREADS, cpus+1) as int
                int queue = property(NXF_QUEUE_SIZE, 1000) as int
                return CustomThreadPool.boundedPool(size, queue)

            case 'unbound':
                int cpus = Runtime.runtime.availableProcessors()
                int size = property(NXF_MAX_THREADS, cpus+1) as int
                return CustomThreadPool.unboundedPool(size)

            default:
                throw new IllegalAccessException("Unknown thread pool type: `$type`")
        }

    }

    @Override
    Pool createPool(boolean daemon) {
        throw new UnsupportedOperationException()
    }

    @Override
    Pool createPool(int numberOfThreads) {
        throw new UnsupportedOperationException()
    }

    @Override
    Pool createPool(boolean daemon, int numberOfThreads) {
        throw new UnsupportedOperationException()
    }

    def property( String name, defValue ) {
        def result = System.getProperty(name)
        return result ?: defValue
    }

}
