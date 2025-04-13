/*
 * Copyright 2013-2024, Seqera Labs
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
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor

import java.nio.file.Path

/**
 * Defines the defaults method for application flow observer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface TraceObserver {

    /**
     * This method is invoked when the flow is going to start
     */
    default void onFlowCreate(Session session) {}

    /**
     * This method is invoked when the flow is going to start
     */
    default void onFlowBegin() {}

    /**
     * This method is invoked when the flow is going to complete
     */
    default void onFlowComplete() {}

    /*
     * Invoked when the process is created.
     */
    default void onProcessCreate( TaskProcessor process ){}

    /*
      * Invoked when all tak have been executed and process ends.
      */
    default void onProcessTerminate( TaskProcessor process ){}

    /**
     * This method when a new task is created and submitted in the nextflow
     * internal queue of pending task to be scheduled to the underlying
     * execution backend
     *
     * @param handler
     * @param trace
     */
    default void onProcessPending(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked before a process run is going to be submitted
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} for the current task.
     */
    default void onProcessSubmit(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run is going to start
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} for the current task.
     */
    default void onProcessStart(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run completes
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} for the current task.
     */
    default void onProcessComplete(TaskHandler handler, TraceRecord trace){}

    /**
     * method invoked when a task execution is skipped because the result is cached (already computed)
     * or stored (due to the usage of `storeDir` directive)
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task
     * @param trace
     *      The trace record for the cached trace. When this event is invoked for a store task
     *      the {@code trace} record is expected to be {@code null}
     */
    default void onProcessCached(TaskHandler handler, TraceRecord trace){}

    /**
     * @return {@code true} whenever this observer requires to collect task execution metrics
     */
    default boolean enableMetrics(){ false }

    /**
     * Method that is invoked, when a workflow fails.
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} for the current task.
     */
    default void onFlowError(TaskHandler handler, TraceRecord trace){}

    /**
     * Method that is invoked when a value is published from a channel.
     *
     * @param value
     */
    default void onWorkflowPublish(Object value){}

    /**
     * Method that is invoke when an output file is published
     * into a `publishDir` folder.
     *
     * @param destination
     *      The destination path at `publishDir` folder.
     */
    default void onFilePublish(Path destination){}

    /**
     * Method that is invoke when an output file is published
     * into a `publishDir` folder.
     *
     * @param destination
     *      The destination path at `publishDir` folder.
     * @param source
     *      The source path at `workDir` folder.
     */
    default void onFilePublish(Path destination, Path source){
        onFilePublish(destination)
    }
}
