/*
 * Copyright 2013-2026, Seqera Labs
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

import java.util.concurrent.ExecutorService
import java.util.function.Supplier

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
/**
 * A configurable factory of the dataflow thread pool
 *
 * See http://niklasschlimm.blogspot.com.es/2012/03/threading-stories-about-robust-thread.html
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CustomPoolFactory implements Supplier<ExecutorService> {

    public final String NXF_POOL_TYPE = 'nxf.pool.type'

    public final String NXF_MAX_THREADS = 'nxf.pool.maxThreads'

    public final String NXF_QUEUE_SIZE = 'nxf.pool.queueSize'

    @Memoized
    @Override
    synchronized ExecutorService get() {

        final defType = Threads.useVirtual() ? 'virtual' : 'default'
        final type = property(NXF_POOL_TYPE, defType)
        switch (type) {
            case 'default':
                int poolSize = Runtime.runtime.availableProcessors() +1
                int maxThreads = property(NXF_MAX_THREADS, 1000) as int
                return CustomThreadPool.defaultPool(poolSize, maxThreads)

            case 'legacy':
                return CustomThreadPool.legacyPool()

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

            case 'virtual':
                log.debug "Creating virtual thread pool"
                return CustomThreadPool.virtualPool()

            default:
                throw new IllegalAccessException("Unknown thread pool type: `$type`")
        }

    }

    def property( String name, defValue ) {
        def result = System.getProperty(name)
        return result ?: defValue
    }

}
