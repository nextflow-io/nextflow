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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.executor.ExecutorConfig
import nextflow.util.Duration

/**
 * Task monitor for the {@code agent} executor. An agent task is an ORCHESTRATOR: its body blocks
 * waiting for the tool sub-tasks it dispatches, which run on the normal (throttled) compute
 * executor. If agents were throttled by the same cpu/capacity budget as those sub-tasks, a set of
 * blocked agents would hold every slot their own sub-tasks need and the run would deadlock.
 *
 * So this monitor does NOT resource-throttle: no cpu/memory accounting (unlike
 * {@link LocalPollingMonitor}) and no queue-capacity cap. Agents are admitted as soon as they are
 * ready; real compute stays bounded by its own executor's monitor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AgentPollingMonitor extends TaskPollingMonitor {

    protected AgentPollingMonitor(Map params) {
        super(params)
    }

    static AgentPollingMonitor create(Session session, ExecutorConfig config, String name) {
        assert session
        assert config
        assert name
        final pollInterval = config.getPollInterval(name, Duration.of('100ms'))
        final dumpInterval = config.getMonitorDumpInterval(name)
        log.debug "Creating agent task monitor for executor '$name' > unbounded (no cpu/capacity throttle); pollInterval: $pollInterval"
        // capacity omitted => 0 => no queue-capacity cap (see TaskPollingMonitor.canSubmit)
        new AgentPollingMonitor(name: name, session: session, config: config, pollInterval: pollInterval, dumpInterval: dumpInterval)
    }

    /**
     * Admit an agent task as soon as it is ready — no cpu/capacity throttling. maxForks
     * ({@code canForkProcess}) is still honoured so an explicit user cap on the agent still applies.
     */
    @Override
    protected boolean canSubmit(TaskHandler handler) {
        handler.canForkProcess() && handler.isReady()
    }
}
