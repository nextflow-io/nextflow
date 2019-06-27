/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.processor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.concurrent.TimeUnit
import com.google.common.util.concurrent.RateLimiter

import nextflow.wr.client.WrRestApi
import nextflow.wr.executor.WrTaskHandler

import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskFault
import nextflow.processor.BatchContext
import nextflow.processor.BatchHandler
import nextflow.util.Throttle
import nextflow.util.Duration
import nextflow.Session
import nextflow.executor.BatchCleanup

/**
 * Monitors the queued tasks, adding them to wr in batches and waiting for their
 * termination
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TaskPollingMonitor by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WrMonitor extends TaskPollingMonitor {

    final private WrRestApi client

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>session: The current {@code Session}
     * <li>client: An instance of WrRestApi (required)
     *
     * @param params
     */
    protected WrMonitor( Map params ) {
        super(params)
        assert params
        assert params.client instanceof WrRestApi
        this.client = params.client as WrRestApi
    }

    static WrMonitor create( Session session, WrRestApi client ) {
        assert session
        assert client

        final String name = 'wr'
        final defPollInterval = Duration.of('1s')
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)
        final int size = session.getQueueSize(name, 1000000)

        log.debug "Creating task monitor for executor '$name' > dumpInterval: $dumpInterval"
        new WrMonitor(
            client: client,
            name: name,
            session: session,
            capacity: size,
            pollInterval: pollInterval,
            dumpInterval: dumpInterval,
        )
    }

    /**
     * Submits the specified tasks for execution adding them to the queue of scheduled tasks
     *
     * @param handler
     *      A List of {@link WrTaskHandler} instance representing the task to be submitted for execution
     */
    protected void submit(WrTaskHandler handler, String id) {
        // we have already done a batched submission by this point, so we're
        // actually just updating state
        handler.submitted(id)
        // note: add the 'handler' into the polling queue *after* the submit operation,
        // this guarantees that in the queue are only jobs successfully submitted
        runningQueue.add(handler)
        // notify task submission
        session.notifyTaskSubmit(handler)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean canSubmit(TaskHandler handler) {
        return true
    }

    /**
     * Loop over the queue of pending tasks and submit all
     *
     * @return The number of tasks submitted for execution
     */
    @Override
    protected int submitPendingTasks() {
        int count = 0
        def itr = getPendingQueue().iterator()
        List<List> toSubmit = []
        Map<String,WrTaskHandler> handlers = [:]
        while( itr.hasNext() ) {
            final handler = itr.next() as WrTaskHandler
            try {
                if( session.isSuccess() ) {
                    itr.remove(); count++   // <-- remove the task in all cases
                    List args = handler.submitArgs()
                    toSubmit << args
                    handlers[args[0] as String] = handler
                }
                else
                    break
            }
            catch ( Throwable e ) {
                handleException(handler, e)
                session.notifyTaskComplete(handler)
            }
        }

        List<Map> jobs = client.add(toSubmit)

        for ( job in jobs ) {
            WrTaskHandler handler = handlers[job."Cmd" as String]
            if (handler) {
                submit(handler, job."Key" as String)
            }
        }

        log.debug("submitPendingTasks looped through $count pending tasks")
        return count
    }

}
