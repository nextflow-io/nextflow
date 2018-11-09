/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
    void onFlowStart(Session session){}

    /**
     * This method is invoked when the flow is going to complete
     */
    void onFlowComplete(){}

    /*
     * Invoked when the process is created.
     */
    void onProcessCreate( TaskProcessor process ){}

    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    void onProcessSubmit(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    void onProcessStart(TaskHandler handler, TraceRecord trace){}

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    void onProcessComplete(TaskHandler handler, TraceRecord trace){}

    /**
     * method invoked when a task execution is skipped because a cached result is found
     * @param handler
     */
    void onProcessCached(TaskHandler handler, TraceRecord trace){}

    /**
     * @return {@code true} whenever this observer requires to collect task execution metrics
     */
    boolean enableMetrics(){ false }

    /**
     * Method that is invoked, when a workflow fails.
     * @param handler
     * @param trace
     */
    void onFlowError(TaskHandler handler, TraceRecord trace){}
}