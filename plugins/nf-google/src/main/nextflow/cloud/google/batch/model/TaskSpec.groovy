/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 * Model a Bath Task Spec
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskSpec {
    List<TaskRunnable> runnables = []
    List<TaskVolume> volumes = []
    ComputeResource computeResource
    /*
     * Maximum duration the task should run
     */
    String maxRunDuration

    /**
     * Maximum number of retries on failures
     * The default, 0, means never retry
     */
    Integer maxRetryCount

    /**
     * Environment variables to set before running the Task.
     * You can set up to 100 environments.
     */
    Map<String,String> environment

    TaskSpec addRunnable(TaskRunnable it) {
        runnables.add(it)
        return this
    }

    TaskSpec withVolumes(List<TaskVolume> volumes) {
        this.volumes = new ArrayList<>(volumes)
        return this
    }

    TaskSpec withComputeResources(ComputeResource it) {
        this.computeResource = it
        return this
    }

    TaskSpec withEnvironment(Map<String,String> env) {
        this.environment = env
        return this
    }


    TaskSpec withMaxRunDuration(Duration duration) {
        this.maxRunDuration = duration ? "${duration.toSeconds()}s" : null
        return this
    }
}
