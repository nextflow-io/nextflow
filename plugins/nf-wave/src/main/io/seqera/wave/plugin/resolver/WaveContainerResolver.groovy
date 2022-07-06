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
import java.util.concurrent.Callable

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import groovy.transform.CompileStatic
import io.seqera.wave.plugin.SubmitContainerTokenResponse
import io.seqera.wave.plugin.WaveClient
import nextflow.Global
import nextflow.Session
import nextflow.container.resolver.ContainerResolver
import nextflow.container.resolver.DefaultContainerResolver
import nextflow.plugin.Priority
import nextflow.processor.TaskRun
/**
 * Implement Wave container resolver logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Priority(-10)  // <-- lower is higher, this is needed to override default provider behavior
class WaveContainerResolver implements ContainerResolver {

    private ContainerResolver defaultResolver = new DefaultContainerResolver()
    private List<String> DOCKER_LIKE = ['docker','podman']
    private final String DOCKER_PREFIX = 'docker://'
    private WaveClient client0


    private Cache<String, SubmitContainerTokenResponse> cache = CacheBuilder<String, SubmitContainerTokenResponse>
                    .newBuilder()
                    .build()


    synchronized protected WaveClient client() {
        if( client0 )
            return client0
        return client0 = new WaveClient( Global.session as Session )
    }

    @Override
    String resolveImage(TaskRun task, String imageName) {
        if( !client().enabled() )
            return defaultResolver.resolveImage(task, imageName)

        final engine = task.processor.executor.isContainerNative()
                ? 'docker'  // <-- container native executor such as AWS Batch are implicitly docker based
                : task.getContainerConfig().getEngine()
        if( engine in DOCKER_LIKE ) {
            final targetImage = defaultResolver.resolveImage(task, imageName)
            return waveContainer(task, targetImage)
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
            final targetImage = waveContainer(task, imageName)
            // then adapt it to singularity format
            return defaultResolver.resolveImage(task,targetImage)
        }
        else {
            // other engine are not supported by wave
            return defaultResolver.resolveImage(task, imageName)
        }
    }

    synchronized String waveContainer(TaskRun task, String container) {
        final bundle = task.getModuleBundle()
        final key = bundle.fingerprint()
        final result = cache.get(key, { client().sendRequest(container, bundle) } as Callable )
        return result.targetImage
    }
}
