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
import nextflow.container.ContainerConfig
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

    private DefaultContainerResolver defaultResolver = new DefaultContainerResolver()
    static final private List<String> DOCKER_LIKE = ['docker','podman','sarus']
    static final private List<String> SINGULARITY_LIKE = ['singularity','apptainer']
    static final private String DOCKER_PREFIX = 'docker://'
    private WaveClient client0

    synchronized protected WaveClient client() {
        if( client0 )
            return client0
        return client0 = new WaveClient( Global.session as Session )
    }

    private String getContainerEngine0(ContainerConfig config) {
        final result = config.isEnabled() ? config.getEngine() : 'docker'
        if( result )
            return result
        // fallback to docker by default
        log.warn "Missing engine in container config - offending value: $config"
        return 'docker'
    }

    @Override
    ContainerInfo resolveImage(TaskRun task, String imageName) {
        if( !client().enabled() )
            return defaultResolver.resolveImage(task, imageName)

        final freeze = client().config().freezeMode()
        final config = task.getContainerConfig()
        final engine = getContainerEngine0(config)
        final singularitySpec = freeze && engine in SINGULARITY_LIKE && !config.canRunOciImage()

        if( engine in DOCKER_LIKE ) {
            // find out the configured image name applying the default resolver
            if( imageName )
                imageName = defaultResolver.resolveImage(task, imageName).getTarget()
            // fetch the wave container image name
            return waveContainer(task, imageName, false)
        }
        else if( engine in SINGULARITY_LIKE ) {
            // remove any `docker://` prefix if any
            if( imageName && imageName.startsWith(DOCKER_PREFIX) )
                imageName = imageName.substring(DOCKER_PREFIX.length())
            // singularity file image use the default resolver
            else if( imageName && (imageName.startsWith('/') || imageName.startsWith('file://') || Files.exists(Path.of(imageName)))) {
                return defaultResolver.resolveImage(task, imageName)
            }
            // fetch the wave container image name
            final image = waveContainer(task, imageName, singularitySpec)
            // when wave returns no info, just default to standard behaviour
            if( !image ) {
                return defaultResolver.resolveImage(task, imageName)
            }
            // oras prefixed container are served directly
            if( image.target.startsWith("oras://") )
                return image
            // otherwise adapt it to singularity format using the target containerInfo to avoid the cache invalidation
            return defaultResolver.resolveImage(task, image.target, image.hashKey)
        }
        else
            throw new IllegalArgumentException("Wave does not support '$engine' container engine")
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
     *      provides a Dockerfile or a Conda recipe or a Spack recipe
     * @return
     *      The container image name returned by the Wave backend or {@code null}
     *      when the task does not request any container or dockerfile to build
     */
    protected ContainerInfo waveContainer(TaskRun task, String container, boolean singularity) {
        final assets = client().resolveAssets(task, container, singularity)
        if( assets ) {
            return client().fetchContainerImage(assets)
        }
        // no container and no dockerfile, wave cannot do anything
        log.trace "No container image or build recipe defined for task ${task.processor.name}"
        return null
    }

    @Override
    boolean isContainerReady(String key) {
        final c=client()
        return c.enabled()
            ? c.isContainerReady(key)
            : defaultResolver.isContainerReady(key)
    }
}
