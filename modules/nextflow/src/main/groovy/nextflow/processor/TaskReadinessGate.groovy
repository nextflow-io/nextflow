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

package nextflow.processor

import groovy.transform.CompileStatic
import org.pf4j.ExtensionPoint

/**
 * Plugin extension point that defers task submission until an external precondition is met.
 *
 * <p>Implementations are invoked once per task by the scheduler, on a managed background
 * thread (virtual thread when available). The task is admitted for submission when this
 * method returns; throw to mark the task as permanently failed and route the cause through
 * the task's {@code errorStrategy}.
 *
 * <p>Implementations may block freely — {@code Thread.sleep}, network calls, long polling.
 * The scheduler thread is never blocked by this call.
 *
 * <p>Implementations <b>must</b> honor {@code Thread.interrupt()} so that task eviction,
 * workflow abort, and the {@code executor.gateMaxWait} backstop can unblock {@code prepare}
 * promptly. Use interruptible primitives ({@code Thread.sleep}, blocking I/O on NIO
 * channels, {@code Future.get}) and propagate {@code InterruptedException}.
 *
 * <p>When multiple gates are registered, a task is admitted only when every gate's
 * {@code prepare} method has returned successfully. Evaluation order across gates is
 * unspecified; all gates start in parallel on the managed executor.
 *
 * <p>Per-process opt-out is available via the {@code hints} directive — gates that wish
 * to support it should check a namespaced hint key (e.g. {@code 'glacier/skip': true}) and
 * return immediately when set. No core mechanism is required.
 */
@CompileStatic
interface TaskReadinessGate extends ExtensionPoint {
    void prepare(TaskHandler handler) throws InterruptedException
}
