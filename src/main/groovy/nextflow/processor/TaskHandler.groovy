/*
 * Copyright (c) 2013, the authors.
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
 * Actions to handle the underlying job running the user task.
 *
 * <p>
 * Note this model the job in the execution facility (i.e. grid, cloud, etc)
 * NOT the *logical* user task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public abstract class TaskHandler {

    enum Status { NEW, SUBMITTED, RUNNING, COMPLETED }

    protected TaskHandler(TaskRun task, TaskConfig taskConfig) {
        this.task = task
        this.taskConfig = taskConfig
    }

    /**
     * The task managed by this handler
     */
    final TaskRun task

    /**
     * The configuration object defined by this task
     */
    final TaskConfig taskConfig

    /**
     * The task managed by this handler
     */
    final TaskRun getTask() { task }

    /**
     * Task current status
     */
    Status status = Status.NEW

    /**
     * Model the start transition from {@code #SUBMITTED} to {@code STARTED}
     */
    abstract boolean checkIfRunning()

    /**
     *  Model the start transition from {@code #STARTED} to {@code TERMINATED}
     */
    abstract boolean checkIfCompleted()

    /**
     * Force the submitted job to quit
     */
    abstract void kill()

    /**
     * Submit the task for execution.
     *
     * Note: the underlying execution platform may schedule it in its own queue
     */
    abstract void submit()


    def void setStatus( Status status ) {

        if ( this.status == status || status == null ) { return }

        this.status = status
    }

    boolean isNew() { return status == Status.NEW }

    boolean isSubmitted() { return status == Status.SUBMITTED }

    boolean isRunning() { return status == Status.RUNNING }

    boolean isCompleted()  { return status == Status.COMPLETED  }

}
