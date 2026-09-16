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
import nextflow.plugin.Plugins
import nextflow.processor.TaskRun

/**
 * Chooses the spec for a task: plugin factories first, then NXF_TASK_HASH_VER, then
 * the current default.
 */
@CompileStatic
class TaskHashSpecResolver {

    static TaskHashSpec resolve(TaskRun task) {
        return resolve(task, Plugins.getPriorityExtensions(TaskHashSpecFactory) ?: Collections.<TaskHashSpecFactory>emptyList())
    }

    static TaskHashSpec resolve(TaskRun task, List<TaskHashSpecFactory> factories) {
        for( TaskHashSpecFactory factory : factories ) {
            final spec = factory.create(task)
            if( spec != null ) {
                return spec
            }
        }
        return defaultSpec()
    }

    static TaskHashSpec defaultSpec() {
        final id = SysEnv.get('NXF_TASK_HASH_VER')
        if( !id ) {
            return StdSpecs.DEFAULT
        }
        return StdSpecs.byId(id)
    }
}
