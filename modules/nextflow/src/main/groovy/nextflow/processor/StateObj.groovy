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