/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

    private T state
    private BlockingDeque events = new LinkedBlockingDeque<>()
    private Thread runner
    private Closure errorHandler

    SimpleAgent(T state) {
        if(state == null)
            throw new IllegalArgumentException("Missing state argument")
        this.state = state
        this.runner = Thread.startDaemon(this.&run)
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
        final retrieve = new RetrieveValueClosure<T>(state)
        events.offerFirst(retrieve)
        return retrieve.getResult()
    }

    T getValue() {
        final retrieve = new RetrieveValueClosure<T>(state)
        events.offer(retrieve)
        return retrieve.getResult()
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
                log.debug "Got an interrupeted exception while polling agent event | ${e.message ?: e}"
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
        Object call(final Object... arguments) {
            result = s0 instanceof Cloneable ? s0.clone() : s0
            sync.countDown()
            return null
        }

        T getResult() {
            try {
                sync.await()
                return (T)result
            }
            catch (InterruptedException e) {
                log.warn "Got an interrupted exception while taking agent result | ${e}"
                return (T)s0
            }
        }

    }
}
