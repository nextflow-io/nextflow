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

import java.util.concurrent.BlockingDeque
import java.util.concurrent.CountDownLatch
import java.util.concurrent.LinkedBlockingDeque
import java.util.concurrent.TimeUnit

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
/**
 * Simple agent that allow to modify and access a mutable state
 * using single lock-free thread
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SimpleAgent<T> {

    /**
     * Max time to wait for the agent runner thread to provide the current state
     */
    static private final Duration GET_VALUE_TIMEOUT = Duration.of('1min')

    private T state
    private BlockingDeque events = new LinkedBlockingDeque<>()
    private Thread runner
    private Closure errorHandler

    SimpleAgent(T state) {
        if(state == null)
            throw new IllegalArgumentException("Missing state argument")
        this.state = state
        this.runner = Threads.start("agent-${state.getClass().getSimpleName()}".toString(), this.&run)
    }

    SimpleAgent onError(@ClosureParams(value = SimpleType, options = ['java.lang.Throwable']) Closure handler) {
        this.errorHandler = handler
        return this
    }

    protected T getState() { state }

    /**
     * Use to modify the state object
     *
     * @param action A closure modifying the agent state
     */
    void send(Closure action) {
        events.offer(action)
    }

    /**
     * Retrieve the current agent state.
     *
     * @return
     *  If the state object implements {@link Cloneable} interface
     *  the cloned state otherwise the state object itself.
     */
    T getQuickValue() {
        if( Thread.currentThread()==runner )
            return currentValue0()
        final retrieve = new RetrieveValueClosure<T>(state)
        events.offerFirst(retrieve)
        return retrieve.getResult()
    }

    T getValue() {
        if( Thread.currentThread()==runner )
            return currentValue0()
        final retrieve = new RetrieveValueClosure<T>(state)
        events.offer(retrieve)
        return retrieve.getResult()
    }

    /**
     * Retrieve the state directly, without going through the events queue. It's meant to be
     * used only when the invoking thread is the agent runner itself, that otherwise would
     * deadlock awaiting for an event that only it can serve.
     */
    @CompileDynamic
    private T currentValue0() {
        return (T)(state instanceof Cloneable ? state.clone() : state)
    }

    protected void run() {
        while(true) {
            try {
                final ev = events.poll(200, TimeUnit.MILLISECONDS)
                if( ev == null )
                    continue

                if( ev instanceof Closure )
                    ev.call()
                else
                    throw new IllegalArgumentException("Invalid agent event object: $ev [${ev.getClass().getName()}]")
            }
            catch (InterruptedException e) {
                log.debug "Got an interrupted exception while polling agent event | ${e.message ?: e}"
                break
            }
            catch(Throwable e) {
                log.debug "Unexpected error while polling agent event | ${e.message ?: e}"
                errorHandler?.call(e)
            }
        }
    }

    @CompileStatic
    private static class RetrieveValueClosure<T> extends Closure {

        private Object s0
        private volatile Object result
        private CountDownLatch sync

        private RetrieveValueClosure(Object value) {
            super(null)
            s0 = value
            sync = new CountDownLatch(1)
        }

        @Override
        @CompileDynamic
        Object call(final Object... arguments) {
            result = s0 instanceof Cloneable ? s0.clone() : s0
            sync.countDown()
            return null
        }

        T getResult() {
            try {
                // note: do not await indefinitely, otherwise a stalled runner thread
                // would hang the invoking thread forever -- see issue #7444
                if( !sync.await(GET_VALUE_TIMEOUT.millis, TimeUnit.MILLISECONDS) ) {
                    log.warn "Timed out awaiting the agent result (>$GET_VALUE_TIMEOUT) -- Returning the current state"
                    return (T)s0
                }
                return (T)result
            }
            catch (InterruptedException e) {
                log.warn "Got an interrupted exception while taking agent result | ${e}"
                return (T)s0
            }
        }

    }
}
