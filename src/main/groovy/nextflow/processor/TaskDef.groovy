/*
 * Copyright (c) 2012, the authors.
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

class TaskDef {

    enum Status { PENDING, RUNNING, TERMINATED }

    /**
     * Task unique id
     */
    def id

    /**
     * Task index
     */
    def index

    /**
     * Task current status
     */
    Status status

    /**
     * The exit code returned by executing the task script
     */
    int exitCode

    /**
     * The return std out
     */
    String output

    /**
     * The directory used to run the task
     */
    File workDirectory

    /**
     * Task start time
     */
    long startTime

    /**
     * Task end time
     */
    long terminateTime


    def void setStatus( Status status ) {
        assert status
        if ( this.status == status ) return

        this.status = status
        if ( status == Status.RUNNING ) {
            startTime = System.currentTimeMillis()
        }
        else if ( status == Status.TERMINATED ) {
            terminateTime = System.currentTimeMillis()
        }
    }

}

