/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.processor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Models a task batch, which executes a set of tasks
 * on the same node.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskBatchRun extends TaskRun {

    List<TaskHandler> children

    @Override
    boolean isContainerEnabled() {
        return false
    }

    void finalize() {
        for( TaskHandler handler : children ) {
            final task = handler.task
            task.exitStatus = exitStatus
            task.error = error
            task.stdout = task.workDir.resolve(TaskRun.CMD_OUTFILE)
            task.stderr = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        }
    }

}
