/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.processor

import java.util.concurrent.Callable
import java.util.concurrent.TimeUnit

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.ThrottlingExecutor
/**
 * Polling monitor class submitting job execution in a parallel manner
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ParallelPollingMonitor extends TaskPollingMonitor {

    private ThrottlingExecutor submitter

    private ThrottlingExecutor reaper

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>name: The name of the executor for which the polling monitor is created
     * <li>session: The current {@code Session}
     * <li>capacity: The maximum number of this monitoring queue
     * <li>pollInterval: Determines how often a poll occurs to check for a process termination
     * <li>dumpInterval: Determines how often the executor status is written in the application log file
     *
     * @param params
     */
    ParallelPollingMonitor(ThrottlingExecutor executor, ThrottlingExecutor reaper, Map params) {
        super(params)
        this.submitter = executor
        this.reaper = reaper
    }

    protected RateLimiter createSubmitRateLimit() {
        // disable the rate limiter managed by the super-class
        // it has be managed by the throttling executor
        return null
    }

    final protected void submit0(TaskHandler handler) {
        super.submit(handler)
    }

    @Override
    protected void submit(TaskHandler handler) {
        // execute task submission in a parallel
        // using an thread-pool via the executor service
        final wrapper = (Callable)new ThrottlingExecutor.Recoverable() {
            @Override protected Object invoke() {
                return submit0(handler)
            }

            // when the submission fails it handles the error
            // condition and notify the task completion (as a failure)
            // note: depending the task error strategy it may try
            // to submit a new task instance
            @Override
            protected void onFailure(Throwable e) {
                if( !session.success )
                    return // ignore error when the session has been interrupted 
                handleException(handler, e)
                session.notifyTaskComplete(handler)
            }
        }

        submitter.submit(wrapper)
    }

    @Override
    protected void cleanup() {
        def tasks = submitter.shutdownNow()
        if( tasks ) log.warn "Execution interrupted -- cleaning up execution pool"
        submitter.awaitTermination(5, TimeUnit.MINUTES)
        // -- now cleanup pending task
        super.cleanup()
        // -- finally delete cleanup executor
        reaper.shutdown()
        reaper.awaitTermination(5, TimeUnit.MINUTES)
    }
}
