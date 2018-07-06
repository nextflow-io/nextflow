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

import com.google.common.util.concurrent.RateLimiter
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.TaskRecoverException
import nextflow.util.Duration
import nextflow.util.RateUnit
import nextflow.util.ThrottlingExecutor
/**
 * Polling monitor class submitting job execution
 * is a parallel manner
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ParallelPollingMonitor extends TaskPollingMonitor {

    private ThrottlingExecutor submitter

    static ParallelPollingMonitor create(Session session, String name, Duration defPollInterval) {
        assert session
        assert name

        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        log.debug "Creating parallel monitor for executor '$name' > pollInterval=$pollInterval; dumpInterval=$dumpInterval"
        new ParallelPollingMonitor(name: name, session: session, pollInterval: pollInterval, dumpInterval: dumpInterval)
    }


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
    protected ParallelPollingMonitor(Map params) {
        super(params)
    }

    @Override
    TaskMonitor start() {
        createExecutorService()
        super.start()
    }

    ThrottlingExecutor getExecutor() { submitter }

    @Override
    protected void cleanup() {
        submitter.shutdownNow()
        super.cleanup()
    }

    private void createExecutorService() {
        final qs = session.getQueueSize(name, 5_000)
        final limit = session.getExecConfigProp(name,'submitRateLimit','50/s') as String
        final size = Runtime.runtime.availableProcessors() * 5

        final opts = new ThrottlingExecutor.Options()
                            .retryOn(TaskRecoverException)
                            .onFailure { Throwable t -> session?.abort(t) }
                            .onRateLimitChange { RateUnit rate -> logRateLimitChange(rate) }
                            .withRateLimit(limit)
                            .withQueueSize(qs)
                            .withPoolSize(size)
                            .withKeepAlive(Duration.of('1 min'))
                            .withAutoThrottle(true)
                            .withMaxRetries(10)
                            .withOptions( getConfigOpts() )

        submitter = ThrottlingExecutor.create(opts)
    }

    @CompileDynamic
    protected Map getConfigOpts() {
        session.config?.executor?.submitter as Map
    }

    /* disable default implementation */
    protected RateLimiter createSubmitRateLimit() { return null }

    protected void logRateLimitChange(RateUnit rate) {
        log.debug "New submission rate limit: $rate"
    }

    @Override
    protected void submit(TaskHandler handler) {
        // execute task submission in a parallel
        // using an thread-pool via the executor service
        submitter.submit({

            try {
                super.submit(handler)
            }
            catch (TaskRecoverException e) {
                throw  e
            }
            catch (Throwable e) {
                log.debug "Handle exception: $e"
                handleException(handler, e)
                session.notifyTaskComplete(handler)
            }

        } as Runnable)
    }
}
