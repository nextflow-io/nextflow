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

package nextflow.trace

import com.google.common.annotations.Beta
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.TaskProcessor
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.trace.event.WorkflowOutputEvent

/**
 * Interface for consuming workflow lifecycle events.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Beta
@CompileStatic
interface TraceObserverV2 {

    /**
     * Invoked when the run is created (before script execution).
     *
     * @param session
     */
    default void onFlowCreate(Session session) {}

    /**
     * Invoked when the workflow DAG is ignited (after script execution).
     */
    default void onFlowBegin() {}

    /**
     * Invoked when the workflow is complete:
     *
     * - all tasks have completed (or aborted)
     * - all files have been published (or aborted)
     * - all workflow onComplete handlers have been executed
     */
    default void onFlowComplete() {}

    /**
     * Invoked when a process is created in the workflow DAG.
     *
     * @param process
     */
    default void onProcessCreate(TaskProcessor process) {}

    /**
     * Invoked when a process is terminated (all process tasks have completed).
     *
     * @param process
     */
    default void onProcessTerminate(TaskProcessor process) {}

    /**
     * Invoked when a task is created and submitted to an executor
     * (but not the execution backend such as the underlying scheduler).
     *
     * @param event
     */
    default void onTaskPending(TaskEvent event) {}

    /**
     * Invoked when a task is submitted by an executor to the
     * underlying execution backend.
     *
     * @param event
     */
    default void onTaskSubmit(TaskEvent event) {}

    /**
     * Invoked when a task is running in the underlying execution backend.
     *
     * @param event
     */
    default void onTaskStart(TaskEvent event) {}

    /**
     * Invoked when a task completes.
     *
     * @param event
     */
    default void onTaskComplete(TaskEvent event) {}

    /**
     * Invoked when a task execution is skipped because the result is cached (already computed)
     * or stored (using the `storeDir` directive).
     *
     * @param event
     */
    default void onTaskCached(TaskEvent event) {}

    /**
     * Whether this observer requires task execution metrics to be collected.
     */
    default boolean enableMetrics() {
        return false
    }

    /**
     * Invoked when a workflow run fails due to a task failure.
     *
     * @param event
     */
    default void onFlowError(TaskEvent event) {}

    /**
     * Invoked when a workflow output is completed.
     *
     * For a given workflow output, this event is guaranteed
     * to be emitted after all files for that output have been
     * published via onFilePublish().
     *
     * @param event
     */
    default void onWorkflowOutput(WorkflowOutputEvent event) {}

    /**
     * Invoked when a file is published by a workflow output
     * or `publishDir` directive.
     *
     * @param event
     */
    default void onFilePublish(FilePublishEvent event) {}

}
