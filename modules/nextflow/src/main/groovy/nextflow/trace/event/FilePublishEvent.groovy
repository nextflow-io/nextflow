/*
 * Copyright 2013-2024, Seqera Labs
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
 * Models a file publish event, which is emitted for each file
 * that is published in a workflow output or publishDir.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Canonical
@CompileStatic
class FilePublishEvent {
    /**
     * The source path.
     */
    Path source
    /**
     * The target path.
     */
    Path target
    /**
     * Labels associated with the published file.
     */
    List<String> labels
}
