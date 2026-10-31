// GPars - Groovy Parallel Systems
//
// Copyright © 2008-11  The original author or authors
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
// Modified by Seqera Labs for Nextflow (2026): a channel read callback instead of an actor message stream.

package groovyx.gpars.dataflow;

import groovy.lang.Closure;
import groovyx.gpars.dataflow.impl.ReadCallback;
import groovyx.gpars.group.PGroup;

/**
 * A helper class enabling the 'whenBound()' functionality of a DataflowVariable and DataflowQueue.
 * Once the channel hands over a value, it runs the supplied closure with the value as a parameter
 * on the thread pool of the parallel group.
 *
 * @author Vaclav Pech, Alex Tkachman
 *         Date: Sep 13, 2009
 */
public final class DataCallback implements ReadCallback {
    private final Closure code;
    private final PGroup parallelGroup;

    /**
     * @param code   The closure to run
     * @param pGroup pGroup The parallel group to join
     */
    public DataCallback(final Closure code, final PGroup pGroup) {
        if (pGroup == null)
            throw new IllegalArgumentException("Cannot create a DataCallback without a parallelGroup parameter");
        if (code.getMaximumNumberOfParameters() > 1)
            throw new IllegalArgumentException("The supplied closure expects more than one argument.");
        this.parallelGroup = pGroup;
        this.code = code;
    }

    /**
     * Schedules the closure with the thread pool.
     * Registers its parallel group with DataflowExpressions for nested 'whenBound' handlers to use the same group.
     */
    @Override
    public void onItem(final Object item) {
        parallelGroup.getThreadPool().execute(new Runnable() {
            @Override
            public void run() {
                Dataflow.activeParallelGroup.set(parallelGroup);
                try {
                    if (code.getMaximumNumberOfParameters() == 1) code.call(item);
                    else code.call();
                } finally {
                    Dataflow.activeParallelGroup.remove();
                }
            }
        });
    }
}
