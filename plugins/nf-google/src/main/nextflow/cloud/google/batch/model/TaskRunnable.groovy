/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
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
