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

package nextflow.trace
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor

/**
 * Defines the defaults method for application flow observer
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
trait TraceObserver {

    /**
     * The is method is invoked when the flow is going to start
     */
    void onFlowCreate(Session session) {}

    /**
     * The is method is invoked when the flow is going to start
     */
    void onFlowBegin() {}

    /**
     * This method is invoked when the flow is going to complete
     */
    void onFlowComplete() {}

    /*
     * Invoked when the process is created.
     */
    void onProcessCreate( TaskProcessor process ){}

    /*
      * Invoked when all tak have been executed and process ends.
      */
    void onProcessTerminate( TaskProcessor process ){}

    /**
     * This method when a new task is created and submitted in the nextflow
     * internal queue of pending task to be scheduled to the underlying
     * execution backend
     *
     * @param handler
     * @param trace
     */
    void onProcessPending(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked before a process run is going to be submitted
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} fot the current task.
     */
    void onProcessSubmit(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run is going to start
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} fot the current task.
     */
    void onProcessStart(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run completes
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} fot the current task.
     */
    void onProcessComplete(TaskHandler handler, TraceRecord trace){}

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
    void onProcessCached(TaskHandler handler, TraceRecord trace){}

    /**
     * @return {@code true} whenever this observer requires to collect task execution metrics
     */
    boolean enableMetrics(){ false }

    /**
     * Method that is invoked, when a workflow fails.
     *
     * @param handler
     *      The {@link TaskHandler} instance for the current task.
     * @param trace
     *      The associated {@link TraceRecord} fot the current task.
     */
    void onFlowError(TaskHandler handler, TraceRecord trace){}
}