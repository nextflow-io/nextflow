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
import nextflow.SysEnv
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHasher
import nextflow.processor.TaskHasherFactory
import nextflow.processor.TaskRun

/**
 * Supplies a {@link VersionedTaskHasher} when the run asked for a task hash version, and
 * abstains otherwise.
 *
 * Abstaining is the whole point: with no version requested this returns {@code null},
 * {@code TaskProcessor.createTaskHasher} falls through to the stock {@link TaskHasher},
 * and the default hashing path is untouched — no spec, no registry, and no risk of
 * moving anyone's cache keys.
 */
@Slf4j
@CompileStatic
class VersionedTaskHasherFactory implements TaskHasherFactory {

    static final String ENV_VAR = 'NXF_TASK_HASH_VER'

    /**
     * @return the requested spec, or {@code null} when none was requested. An unknown id
     *      is an error rather than a silent fallback: a run hashed under a version other
     *      than the one it asked for produces cache keys nobody can account for.
     */
    static TaskHashSpec requestedSpec() {
        final id = SysEnv.get(ENV_VAR)
        if( !id ) {
            return null
        }
        return StdSpecs.byId(id)
    }

    @Override
    TaskHasher create(TaskRun task) {
        final spec = requestedSpec()
        if( spec == null ) {
            return null
        }
        log.debug "Task: ${task.lazyName()} > hashing under spec ${spec.id} (${spec.fingerprint()})"
        return new VersionedTaskHasher(task, spec)
    }
}
