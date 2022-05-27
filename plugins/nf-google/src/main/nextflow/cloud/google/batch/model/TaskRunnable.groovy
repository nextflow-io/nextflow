/*
 * Copyright 2022, Google Inc.
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

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model a Google Batch task configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class TaskRunnable {
    TaskContainer container

    // Normally, a non−zero exit status causes the Task // execution of other Runnables to continue instead.
    Boolean ignoreExitStatus

    /**
     * This flag allows a Runnable to continue running in the background while the
     * Task executes subsequent Runnables. This is useful to provide services to
     * other Runnables (or to provide debugging support tools like SSH servers).
     */
    Boolean background

    /**
     * By default, after a Runnable fails, no further Runnable are executed. This
     * flag indicates that this Runnable must be run even if the Task has already
     * failed. This is useful for Runnables that copy output files off of the VM
     * or for debugging
     *
     * The always_run flag does not override the Task's overall max_run_duration.
     * If the max_run_duration has expired then no further Runnables will execute,
     * not even always_run Runnables.
     */
    Boolean alwaysRun

    TaskRunnable withContainer(TaskContainer it) {
        this.container = it
        return this
    }
}
