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

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
/**
 * Models Seqera Intelligent Compute scheduler metadata for Nextflow execution.
 *
 * <p>The {@code enabled} flag is derived from the run configuration (i.e. whether the
 * {@code seqera} executor is selected) so that it is available at workflow create/begin
 * time, before any executor is registered. It mirrors {@link WaveMetadata} and
 * {@link FusionMetadata}.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
class SchedulerMetadata {

    /**
     * Name of the executor backed by the Seqera Intelligent Compute scheduler.
     * Matches the {@code ServiceName} of the {@code SeqeraExecutor} in the {@code nf-seqera} plugin
     * (referenced as a literal to avoid a core → plugin dependency).
     */
    private static final String SEQERA_EXECUTOR = 'seqera'

    final boolean enabled

    /**
     * The scheduler-assigned run identifier.
     *
     * Volatile because it is written by the {@code SeqeraExecutor} on the executor thread
     * (lazily, on the first task submission) and read by the {@code TowerObserver} reporting
     * thread. It is propagated to Platform via the top-level {@code schedulerRunId} field on
     * trace progress and heartbeat requests — not as part of the serialized workflow object.
     */
    volatile String runId

    SchedulerMetadata(Session session) {
        this( resolveExecutor(session) )
    }

    SchedulerMetadata(String executorName) {
        this.enabled = executorName == SEQERA_EXECUTOR
    }

    /**
     * Render the metadata for serialization on the workflow object of begin/complete requests.
     *
     * <p>{@link #runId} is included only when assigned: it is unset at begin (it is assigned on
     * the first task submission) so the begin request carries just {@code enabled}, while the
     * complete request also carries the {@code runId}. The id is additionally propagated — at the
     * earliest moment it exists — via the top-level {@code schedulerRunId} field on trace progress
     * requests.
     */
    Map toMap() {
        final result = new LinkedHashMap(2)
        result.enabled = enabled
        if( runId != null )
            result.runId = runId
        return result
    }

    /**
     * Resolve the run-level executor name from the configuration, following the same
     * precedence used by {@code ExecutorFactory}: {@code process.executor}, then
     * {@code executor.name}, then the {@code NXF_EXECUTOR} environment variable.
     *
     * <p>Note: per-process executor overrides (e.g. {@code withName:...{ executor='seqera' }})
     * are not considered here; the flag reflects the run-level executor selection.
     */
    private static String resolveExecutor(Session session) {
        final process = session?.config?.process as Map
        if( process?.executor )
            return process.executor.toString()
        final executor = session?.config?.executor as Map
        if( executor?.name )
            return executor.name.toString()
        return SysEnv.get('NXF_EXECUTOR')
    }
}