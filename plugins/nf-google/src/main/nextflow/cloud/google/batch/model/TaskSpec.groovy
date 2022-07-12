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
import nextflow.util.Duration

/**
 * Model a Google Batch Task specification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true, ignoreNulls = true, includePackage = false)
class TaskSpec {

    List<TaskRunnable> runnables = []

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

    // Lifecycle management schema when any task in a task group is failed.
    // The valid size of lifecycle policies are [0, 10].
    // For each lifecycle policy, when the condition is met,
    // the action in that policy will be executed.
    // If there are multiple policies that the task execution result matches,
    // we use the action from the first matched policy. If task execution result
    // does not meet with any of the defined lifecycle policy, we consider it as
    // the default policy. Default policy means if the exit code is 0, exit task.
    // If task ends with non-zero exit code, retry the task with max_retry_count.
    List<LifecyclePolicy> lifecyclePolicies;

    /**
     * Environment variables to set before running the Task.
     * You can set up to 100 environments.
     */
    Map<String,String> environment

    /*
     * Volumes to mount before running Tasks using this TaskSpec.
     */
    List<TaskVolume> volumes = []

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
