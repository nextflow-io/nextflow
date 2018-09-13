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