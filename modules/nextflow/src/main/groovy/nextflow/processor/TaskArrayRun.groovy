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
import nextflow.container.ContainerConfig
import nextflow.executor.TaskArrayExecutor

/**
 * Models a task array, which submits a collection of independent
 * tasks with a single submit script.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskArrayRun extends TaskRun {

    List<TaskHandler> children

    int getArraySize() {
        children.size()
    }

    @Override
    ContainerConfig getContainerConfig() {
        final config = super.getContainerConfig()
        final envWhitelist = config.getEnvWhitelist() ?: []
        final executor = (TaskArrayExecutor)processor.getExecutor()
        envWhitelist << executor.getArrayIndexName()
        config.put('envWhitelist', envWhitelist)
        return config
    }

    @Override
    boolean isContainerEnabled() {
        return false
    }

    @Override
    final boolean isArray() {
        return true
    }

}
