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
 * Model Batch job notification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class JobNotification {

    static class Message {
        enum Type {
            // Unspecified.
            TYPE_UNSPECIFIED(0),

            // Notify users that the job state has changed.
            JOB_STATE_CHANGED(1),

            // Notify users that the task state has changed.
            TASK_STATE_CHANGED(2);

            private final int value;

            Type(int value) { this.value = value }
        }

    }

    // The Pub/Sub topic where notifications like the job state changes
    // will be published. This topic should be an existing topic in the same
    // project with the job and billings will be charged to this project. If no
    // topic is specified, there will be no Pub/Sub messages sent. Topic
    // format is `projects/{project}/topics/{topic}`.
    String pubsubTopic

    // The message caters the message attributes configuration will to be sent
    // to this Pub/Sub topic. Without this field, there is no message being sent
    // by default.
    Message message

}
