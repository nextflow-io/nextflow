// GPars - Groovy Parallel Systems
//
// Copyright © 2013  The original author or authors
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

// Modified by Seqera Labs for Nextflow (2026): added license header.

package groovyx.gpars.dataflow;

import groovy.lang.Closure;
import groovyx.gpars.actor.impl.MessageStream;
import groovyx.gpars.dataflow.expression.DataflowExpression;

/**
 * Created with IntelliJ IDEA.
 * User: Vaclav
 * Date: 5.7.13
 * Time: 22:17
 * To change this template use File | Settings | File Templates.
 */
public interface SelectableChannel<T> {
    /**
     * Send all pieces of data bound in the future to the provided stream when it becomes available.     *
     *
     * @param closure closure to execute when data becomes available. The closure should take at most one argument.
     */
    <V> void wheneverBound(Closure<V> closure);

    /**
     * Send all pieces of data bound in the future to the provided stream when it becomes available.
     *
     * @param stream stream where to send result
     */
    void wheneverBound(MessageStream stream);

    /**
     * Retrieves the value at the head of the buffer. Returns null, if no value is available.
     *
     * @return The value bound to the DFV at the head of the stream or null
     * @throws InterruptedException If the current thread is interrupted
     */
    @SuppressWarnings({"ClassReferencesSubclass"}) DataflowExpression<T> poll() throws InterruptedException;
}
