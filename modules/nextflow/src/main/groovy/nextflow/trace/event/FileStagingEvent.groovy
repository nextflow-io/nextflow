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
 * Models a file staging event.
 *
 * @author Robrecht Cannoodt <robrecht@data-intuitive.com>
 */
@Canonical
@CompileStatic
class FileStagingEvent {
    /**
     * The original source path (e.g., remote URL or path).
     */
    Path source
    /**
     * The target staged path in the work directory.
     */
    Path target
}
