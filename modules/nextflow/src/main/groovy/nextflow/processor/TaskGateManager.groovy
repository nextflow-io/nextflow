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

import java.util.concurrent.Callable
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessException
import nextflow.plugin.Plugins
import nextflow.util.CustomThreadFactory
import nextflow.util.Threads

/**
 * Manages the lifecycle of {@link TaskReadinessGate} extensions on behalf of a
 * {@link TaskPollingMonitor}: discovers registered gates via PF4J, runs each gate's
 * {@code prepare} call on a managed virtual-thread executor, tracks the resulting
 * futures per task, and surfaces completion to the monitor via {@link #isReady}.
 *
 * <p>The monitor delegates to this class from {@code schedule}, {@code canSubmit},
 * and {@code evict}. When no plugin registers a gate, all operations are no-ops and
 * no executor is created.
 *
 * @author Rob Syme <rob.syme@seqera.io>
 */
@Slf4j
@CompileStatic
class TaskGateManager {

    private final List<TaskReadinessGate> gates

    private ExecutorService executor

    private final ConcurrentMap<TaskHandler, List<Future<?>>> futuresByHandler = new ConcurrentHashMap<>()

    TaskGateManager(Session session) {
        this(session, Plugins.getExtensions(TaskReadinessGate))
    }

    @PackageScope
    TaskGateManager(Session session, List<TaskReadinessGate> gates) {
        this.gates = gates
        if( gates ) {
            this.executor = Threads.useVirtual()
                ? Executors.newVirtualThreadPerTaskExecutor()
                : Executors.newCachedThreadPool(new CustomThreadFactory('TaskReadinessGate'))
            session?.onShutdown { executor.shutdownNow() }
            log.debug "Registered ${gates.size()} task readiness gate(s): ${gates*.class*.simpleName}"
        }
    }

    /**
     * Submit each registered gate's {@code prepare} call for the given handler. The
     * caller is responsible for serializing concurrent calls for the same handler
     * (in {@code TaskPollingMonitor} this happens under {@code pendingLock}).
     */
    void submit(TaskHandler handler) {
        if( !gates ) return

        final futures = new ArrayList<Future<?>>(gates.size())
        for( TaskReadinessGate g : gates ) {
            final gate = g   // capture in a fresh local for the async closure
            futures << executor.submit({ gate.prepare(handler) } as Callable)
        }
        futuresByHandler.put(handler, futures)
    }

    /**
     * Return {@code true} when every gate's future for this handler has completed
     * successfully. Returns {@code false} while at least one future is still running.
     * Throws {@link ProcessException} when any gate has failed; identity-preserves
     * {@code ProcessException} causes and wraps everything else so retry markers
     * carried on the cause reach {@link TaskProcessor#resumeOrDie} intact.
     */
    boolean isReady(TaskHandler handler) {
        final futures = futuresByHandler.get(handler)
        if( !futures ) return true

        for( f in futures ) {
            if( !f.isDone() ) return false
            if( f.isCancelled() ) {
                futures*.cancel(true)
                futuresByHandler.remove(handler)
                throw new ProcessException("Task readiness gate was cancelled for task '${handler.task.name}'")
            }
            try { f.get() }
            catch( ExecutionException e ) {
                // cancel peer gates so their work doesn't outlive the failing task
                futures*.cancel(true)
                futuresByHandler.remove(handler)
                final cause = e.cause
                if( cause instanceof ProcessException )
                    throw (ProcessException) cause
                throw new ProcessException("Task readiness gate failed for task '${handler.task.name}'", cause ?: e)
            }
            catch( InterruptedException e ) {
                // the polling thread itself was interrupted (e.g. session shutdown);
                // restore the flag so subsequent blocking calls fail fast and leave the
                // handler in futuresByHandler so eviction can clean up on the next pass
                Thread.currentThread().interrupt()
                return false
            }
        }
        futuresByHandler.remove(handler)
        return true
    }

    /**
     * Cancel any in-flight gate futures for the given handler. Idempotent — safe to
     * call whether or not gate state exists for the handler.
     */
    void evict(TaskHandler handler) {
        futuresByHandler.remove(handler)?.each { it.cancel(true) }
    }

    /** Visible for testing: snapshot of currently tracked handlers. */
    @PackageScope
    Set<TaskHandler> trackedHandlers() {
        Collections.unmodifiableSet(futuresByHandler.keySet())
    }

    /** Visible for testing: futures for a given handler, or null. */
    @PackageScope
    List<Future<?>> futuresFor(TaskHandler handler) {
        futuresByHandler.get(handler)
    }
}
