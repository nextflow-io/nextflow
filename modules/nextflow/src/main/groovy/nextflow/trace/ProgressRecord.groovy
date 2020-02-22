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
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Holds process execution progress stats
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
class ProgressRecord implements Cloneable {

    final int index
    final String name     // process name
    String hash     // current task hash
    String taskName // current task name
    int pending     // number of new tasks ready to be submitted
    int submitted   // number of tasks submitted for execution not yet started
    int running     // number of tasks whose execution started
    int succeeded   // number of tasks whose execution completed successfully
    int cached      // number of tasks whose execution
    int failed
    int aborted
    int stored
    int ignored
    int retries
    boolean terminated
    boolean errored

    int loadCpus
    long loadMemory

    int peakRunning
    int peakCpus
    long peakMemory

    ProgressRecord() {}

    ProgressRecord(int processId, String processName) {
        this.index = processId
        this.name = processName
        this.taskName = processName
    }

    int getTotalCount() {
        pending+ submitted+ running+
           succeeded+ failed+ cached+ stored
    }

    int getCompletedCount() {
        succeeded+ failed+ cached+ stored
    }

    @Override
    synchronized ProgressRecord clone() {
        return (ProgressRecord)super.clone()
    }

}