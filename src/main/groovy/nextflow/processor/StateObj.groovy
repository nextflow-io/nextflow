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

package nextflow.processor

import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * This object keeps track of a process termination state. This is needed to send out a poison-pill
 * on the output channels, when the process is finished i.e. all tasks have been executed and, at least one,
 * poison-pill has arrived on its input channel.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@EqualsAndHashCode
@PackageScope
class StateObj implements Serializable, Cloneable {

    private int submitted
    private int completed
    private boolean poisoned
    private String name

    StateObj( String name = null ) {
        this.name = name
    }

    /**
     * Increment the number of submitted tasks for execution
     */
    void incSubmitted() {
        if( poisoned )
            log.debug "Oops.. Cannot process more messages after Poison-Pill was received"
        else
            submitted++
    }

    /**
     * Increment the number of completed tasks
     */
    void incCompleted() {
        if( completed >= submitted ) {
            log.debug "Oops.. Processed messages ($submitted) should not overcome received messages ($submitted) count"
        }
        completed++
    }

    void poison() {
        log.trace "<$name> State before poison: $this"
        poisoned = true
    }

    boolean isFinished() {
        poisoned && (submitted == completed)
    }

    StateObj clone() {
        def result = (StateObj)super.clone()
        result.submitted = this.submitted
        result.completed = this.completed
        result.poisoned = this.poisoned
        result.name = this.name
        return result
    }

    String toString() {
        "${getClass().simpleName}[submitted: $submitted; completed: $completed; poisoned: $poisoned ]"
    }

}