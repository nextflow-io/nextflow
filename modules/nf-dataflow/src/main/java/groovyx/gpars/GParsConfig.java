// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2012  The original author or authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Modified by Seqera Labs for Nextflow (2026): removed the timer factory and the shutdown, used by actors, GParsPool and GParsExecutorsPool.

package groovyx.gpars;

import groovyx.gpars.scheduler.Pool;
import groovyx.gpars.scheduler.ResizeablePool;
import groovyx.gpars.util.PoolFactory;

/**
 * Enables to specify custom thread pools to run GPars in hosted environments, such as GAE
 *
 * @author Vaclav Pech
 */
public final class GParsConfig {
    private static volatile PoolFactory poolFactory;
    private static volatile boolean poolFactoryFlag = false;

    public static synchronized void setPoolFactory(final PoolFactory pool) {
        if (poolFactoryFlag)
            throw new IllegalArgumentException("The pool factory cannot be altered at this stage. It has already been set before.");
        poolFactoryFlag = true;
        poolFactory = pool;
    }

    public static PoolFactory getPoolFactory() {
        return poolFactory;
    }

    /**
     * If a pool factory has been set, it will be used to create a new thread pool.
     * Otherwise a new instance of ResizeablePool will be returned.
     *
     * @return A thread pool instance to use for default parallel groups (dataflow, agents)
     */
    public static Pool retrieveDefaultPool() {
        if (poolFactory != null) return poolFactory.createPool();
        return new ResizeablePool(true, 1);
    }
}
