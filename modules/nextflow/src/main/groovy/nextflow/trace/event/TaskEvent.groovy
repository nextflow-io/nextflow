/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.trace.event

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceRecord

/**
 * Models a task lifecycle event.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Canonical
@CompileStatic
class TaskEvent {
    /**
     * The task handler for the given task.
     */
    TaskHandler handler
    /**
     * The trace record for the given task.
     */
    TraceRecord trace
}
