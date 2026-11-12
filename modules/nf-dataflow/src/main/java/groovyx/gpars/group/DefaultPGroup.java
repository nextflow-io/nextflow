// GPars - Groovy Parallel Systems
//
// Copyright © 2008-2011, 2013  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): removed the actor examples from the documentation.

package groovyx.gpars.group;

import groovyx.gpars.GParsConfig;
import groovyx.gpars.scheduler.DefaultPool;
import groovyx.gpars.scheduler.Pool;

/**
 * Provides logical grouping for agents and dataflow operators. Each group has an underlying thread pool,
 * which will perform actions on behalf of the users belonging to the group. Agents and operators created
 * through the {@code agent} or {@code operator} methods of a group automatically belong to it.
 *
 * @author Vaclav Pech
 */
public final class DefaultPGroup extends PGroup {

    /**
     * Creates a group for agents and operators, which will share the supplied thread pool.
     *
     * @param threadPool The thread pool to use for the group
     */
    public DefaultPGroup(final Pool threadPool) {
        super(threadPool);
    }

    /**
     * Creates a group for agents and operators, which will share a common daemon thread pool.
     */
    public DefaultPGroup() {
        super(GParsConfig.getPoolFactory() == null ? new DefaultPool(true) : GParsConfig.getPoolFactory().createPool(true));
    }

    /**
     * Creates a group for agents and operators, which will share a common daemon thread pool.
     *
     * @param poolSize The initial size of the underlying thread pool
     */
    public DefaultPGroup(final int poolSize) {
        super(GParsConfig.getPoolFactory() == null ? new DefaultPool(true, poolSize) : GParsConfig.getPoolFactory().createPool(true, poolSize));
    }
}
