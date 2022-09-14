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

package io.seqera.wave.plugin.resolver

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.wave.plugin.WaveClient
import nextflow.Global
import nextflow.Session
import nextflow.container.resolver.ContainerInfo
import nextflow.container.resolver.ContainerResolver
import nextflow.container.resolver.DefaultContainerResolver
import nextflow.plugin.Priority
import nextflow.processor.TaskRun
/**
 * Implement Wave container resolve logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Priority(-10)  // <-- lower is higher, this is needed to override default provider behavior
class WaveContainerResolver implements ContainerResolver {

    private ContainerResolver defaultResolver = new DefaultContainerResolver()
    private List<String> DOCKER_LIKE = ['docker','podman']
    private final String DOCKER_PREFIX = 'docker://'
    private WaveClient client0

    synchronized protected WaveClient client() {
        if( client0 )
            return client0
        return client0 = new WaveClient( Global.session as Session )
    }

    @Override
    ContainerInfo resolveImage(TaskRun task, String imageName) {
        if( !client().enabled() )
            return defaultResolver.resolveImage(task, imageName)

        if( !imageName ) {
            // when no image name is provider the module bundle should include a
            // Dockerfile or a Conda recipe to build an image on-fly with an
            // automatically assigned name
            return waveContainer(task, null)
        }

        final engine = task.processor.executor.isContainerNative()
                ? 'docker'  // <-- container native executor such as AWS Batch are implicitly docker based
                : task.getContainerConfig().getEngine()
        if( engine in DOCKER_LIKE ) {
            final image = defaultResolver.resolveImage(task, imageName)
            return waveContainer(task, image.target)
        }
        else if( engine=='singularity' ) {
            // remove any `docker://` prefix if any
            if( imageName.startsWith(DOCKER_PREFIX) )
                imageName = imageName.substring(DOCKER_PREFIX.length())
            // singularity file image use the default resolver
            else if( imageName.startsWith('/') || imageName.startsWith('file://') || Files.exists(Path.of(imageName))) {
                return defaultResolver.resolveImage(task, imageName)
            }
            // fetch the wave container name
            final image = waveContainer(task, imageName)
            // then adapt it to singularity format
            return defaultResolver.resolveImage(task, image.target)
        }
        else {
            // other engine are not supported by wave
            return defaultResolver.resolveImage(task, imageName)
        }
    }

    /**
     * Given the target {@link TaskRun} and container image name
     * creates a {@link io.seqera.wave.plugin.WaveAssets} object which holds
     * the corresponding resources to submit container request to the Wave backend
     *
     * @param task
     *      An instance of {@link TaskRun} task representing the current task
     * @param container
     *      The container image name specified by the task. Can be {@code null} if the task
     *      provides a Dockerfile or a Conda recipe
     * @return
     *      The container image name returned by the Wave backend or {@code null}
     *      when the task does not request any container or dockerfile to build
     */
    protected ContainerInfo waveContainer(TaskRun task, String container) {
        final assets = client().resolveAssets(task, container)
        if( assets ) {
            return client().fetchContainerImage(assets)
        }
        // no container and no dockerfile, wave cannot do anything
        log.trace "No container image or build recipe defined for task ${task.processor.name}"
        return null
    }

}
