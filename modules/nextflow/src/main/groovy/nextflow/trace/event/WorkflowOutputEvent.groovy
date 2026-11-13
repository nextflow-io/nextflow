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

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic

/**
 * Models a workflow output event, which is emitted for each
 * workflow output when it is completed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Canonical
@CompileStatic
class WorkflowOutputEvent {
    /**
     * The name of the workflow output.
     */
    String name
    /**
     * The value of the workflow output.
     *
     * If the source is a dataflow channel. this value is an unordered
     * collection of the published values from the channel. If the source
     * is a dataflow value, this value is the published value.
     *
     * If the index file was enabled, this value is null.
     */
    Object value
    /**
     * The optional index file for the workflow output.
     */
    Path index
}
