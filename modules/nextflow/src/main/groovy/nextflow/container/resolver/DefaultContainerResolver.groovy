/*
 * Copyright 2020-2022, Seqera Labs
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
 *
 */

package nextflow.container.resolver

import groovy.transform.CompileStatic
import nextflow.container.ContainerHandler
import nextflow.processor.TaskRun
/**
 * Given a container image name resolves it to target image name.
 *
 * The resolver takes care to prepend the target registry name defined
 * in the nextflow config or adapt the container image name depending
 * the configured container engine
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
final class DefaultContainerResolver implements ContainerResolver {

    @Override
    ContainerInfo resolveImage(TaskRun task, String imageName) {
        if( !imageName ) {
            // no image given, just return null
            return ContainerInfo.EMPTY
        }

        final cfg = task.getContainerConfig()
        final handler = new ContainerHandler(cfg, task.processor.executor)
        final ret = handler.normalizeImageName(imageName)
        return new ContainerInfo(imageName, ret, ret)
    }

}
