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

package nextflow.executor.local

import java.util.concurrent.ExecutorService

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.SupportedScriptTypes
import nextflow.processor.AgentPollingMonitor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType

/**
 * Local executor for {@code agent} tasks. Agents are in-JVM ORCHESTRATORS: their (native) body
 * blocks waiting for the tool sub-tasks it dispatches, which run on the normal {@link LocalExecutor}.
 * Running agents on the standard local executor deadlocks on a small machine because a blocked agent
 * holds a cpu/capacity slot in {@link nextflow.processor.LocalPollingMonitor} that its own sub-task
 * needs.
 *
 * The separation is on TWO axes, and both are required:
 *
 * <ul>
 *   <li><b>Admission</b> — an {@link AgentPollingMonitor} that does not resource-throttle, so a
 *       blocked agent never occupies a cpu slot its tool sub-tasks require.
 *   <li><b>Threads</b> — its own orchestration pool ({@link Session#getAgentExecService}) rather
 *       than the session's execution pool. Sharing one pool puts a dependency edge between two
 *       members of the same pool, so enough blocked agents leave no thread to run the sub-tasks
 *       that would release them. No pool SIZE fixes that, because agents are admitted without a
 *       capacity cap; only a separate pool does.
 * </ul>
 *
 * Real compute (the sub-tasks) stays throttled on the standard local executor, and the dependency
 * crosses the pool boundary in one direction only — orchestration waits on execution, never the
 * reverse — so there is no cycle to close.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@SupportedScriptTypes( [ScriptType.GROOVY] )
class AgentExecutor extends LocalExecutor {

    @Override
    protected TaskMonitor createTaskMonitor() {
        return AgentPollingMonitor.create(session, config, name)
    }

    /**
     * Agents run on the orchestration pool, never on the execution pool they dispatch into.
     */
    @Override
    ExecutorService getExecService() {
        return session.getAgentExecService()
    }

    /**
     * An in-JVM agent runs no wrapper script, so nothing would stage its declared input files.
     * {@link AgentTaskHandler} materializes them into the work dir before the body runs, which is
     * what makes a typed `Path` input mean the same thing here as it does under a container.
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir

        if( task.type == ScriptType.GROOVY )
            return new AgentTaskHandler(task, this)
        return super.createTaskHandler(task)
    }
}
