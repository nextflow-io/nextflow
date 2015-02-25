/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskMonitor {

    /**
     * Put a {@code TaskHandler} instance into the queue of tasks to be processed.
     * Subclasses provide concrete logic to manage them.
     *
     * Note: depending the task scheduling implementation this operation may be blocking
     * to await there's free space in the tasks queue
     *
     * @param handler {@code TaskHandler} instance
     */
    def void put(TaskHandler handler)

    /**
     * Remove the {@code TaskHandler} instance from the queue of tasks ro be processed
     *
     * @param handler A not null {@code TaskHandler} instance
     */
    def boolean drop(TaskHandler handler)

    /**
     * Start the monitoring activity for the queued tasks
     * @return The instance itself, useful to chain methods invocation
     */
    def TaskMonitor start()

    /**
     * Notify when a task terminates
     */
    void signal()
}
