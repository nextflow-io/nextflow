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
package nextflow.processor.hash

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHasher
import nextflow.processor.TaskHasherFactory
import nextflow.processor.TaskRun

/**
 * Supplies a {@link SpecTaskHasher} when the run asked for a task hash version, and
 * abstains otherwise.
 *
 * Abstaining is the whole point: with no version requested this returns {@code null},
 * {@code TaskProcessor.createTaskHasher} falls through to the stock {@link TaskHasher},
 * and the default hashing path is untouched — no spec, no registry, and no risk of
 * moving anyone's cache keys.
 */
@Slf4j
@CompileStatic
class SpecTaskHasherFactory implements TaskHasherFactory {

    @Override
    TaskHasher create(TaskRun task) {
        final spec = TaskHashSpecResolver.requestedSpec()
        if( spec == null ) {
            return null
        }
        log.debug "Task: ${task.lazyName()} > hashing under spec ${spec.id} (${spec.fingerprint()})"
        return new SpecTaskHasher(task, spec)
    }
}
