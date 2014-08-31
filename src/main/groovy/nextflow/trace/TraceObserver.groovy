/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
interface TraceObserver {

    /**
     * The is method is invoked when the flow is going to start
     */
    void onFlowStart(Session session)

    /**
     * This method is invoked when the flow is going to complete
     */
    void onFlowComplete()

    /**
     * TODO This method is invoked when an error is raised
     *
     * @param error
     */
    void onFlowError(Throwable error)

    /*
     * Invoked when the process is created.
     */
    void onProcessCreate( TaskProcessor process )

    /**
     * Invoked when the process has terminated and is going to be destroy
     * @param process
     */
    void onProcessDestroy( TaskProcessor process )

    /**
     * This method is invoked before a process run is going to be submitted
     * @param handler
     */
    void onProcessSubmit(TaskHandler handler)

    /**
     * This method is invoked when a process run is going to start
     * @param handler
     */
    void onProcessStart(TaskHandler handler)

    /**
     * This method is invoked when a process run completes
     * @param handler
     */
    void onProcessComplete(TaskHandler handler)

    /**
     * This method is when a process run raises an error
     * @param handler
     * @param error
     */
    void onProcessError(TaskHandler handler, Throwable error)


}