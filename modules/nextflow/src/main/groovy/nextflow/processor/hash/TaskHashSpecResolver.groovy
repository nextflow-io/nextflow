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

/**
 * Resolves the task hash version a run asked for.
 *
 * Returning {@code null} is the normal case and means "no version requested", which
 * leaves the run on the inherited {@link nextflow.processor.TaskHasher}. An unknown id
 * is an error rather than a silent fallback: a run that asked to be hashed under a
 * particular version and was quietly hashed under another would produce a cache whose
 * keys nobody can account for.
 */
@CompileStatic
class TaskHashSpecResolver {

    static final String ENV_VAR = 'NXF_TASK_HASH_VER'

    /**
     * @return the requested spec, or {@code null} when none was requested.
     * @throws IllegalArgumentException if a version was requested but is not known.
     */
    static TaskHashSpec requestedSpec() {
        final id = SysEnv.get(ENV_VAR)
        if( !id ) {
            return null
        }
        return StdSpecs.byId(id)
    }
}
