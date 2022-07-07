/*
 * Copyright 2020, Seqera Labs
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
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic

/**
 * Model batch job status
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class JobStatus {

    enum State {
        STATE_UNSPECIFIED,

        // Job is submitted into a ResourcePool and waiting
        // for resource allocation.
        QUEUED,

        // Job is scheduled to run as soon as resource allocation is ready.
        // The resource allocation may happen at a later time but with a high
        // chance to succeed.
        SCHEDULED,

        // Resource allocation has been successful. At least one Task in the Job is
        // RUNNING.
        RUNNING,

        // All Tasks in the Job have finished successfully.
        SUCCEEDED,

        // At least one Task in the Job has failed.
        FAILED,

        // The Job will be deleted, but has not been deleted yet. Typically this is
        // because resources used by the Job are still being cleaned up.
        DELETION_IN_PROGRESS
    }

    static class StatusEvent {
        // Type of the event.
        String type

        // Description of the event.
        String description

        // The time this event occurred.
        String eventTime

        // Task Execution
        TaskExecution taskExecution
    }

    static class TaskExecution {
        // When task is completed as the status of FAILED or SUCCEEDED,
        // exit code is for one task execution result, default is 0 as success.
        Integer exitCode
    }

    static class TaskGroupStatus {
        // Count of task in each state in the TaskGroup.
        // The map key is task state name.
        Map<String,Integer> counts
    }

    // Job state
    State state

    // Job status events
    List<StatusEvent> statusEvents

    // Aggregated task status for each TaskGroup in the Job.
    // The map key is TaskGroup ID.
    Map<String, TaskGroupStatus> taskGroups

    // The duration of time the Job is in status
    // RUNNING. Once the Job completes (i.e. the Job status is either
    // SUCCEEDED/FAILED) the run duration represents the time it took the Job
    // to complete.
    //  google.protobuf.Duration run_duration
}
