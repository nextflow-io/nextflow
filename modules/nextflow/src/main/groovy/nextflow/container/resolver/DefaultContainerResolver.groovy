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

import nextflow.container.ContainerHandler
import nextflow.processor.TaskRun

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultContainerResolver implements ContainerResolver {

    @Override
    String resolveImage(TaskRun task, String imageName) {
        final cfg = task.getContainerConfig()
        final handler = new ContainerHandler(cfg, task.processor.executor)
        final result = handler.normalizeImageName(imageName)

        final proxy = System.getenv('NXF_PROXY_REG')
        return proxy ? ContainerHandler.proxyReg(proxy, result) : result
    }
}
