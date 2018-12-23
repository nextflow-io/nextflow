/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
