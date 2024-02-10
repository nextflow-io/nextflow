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

package io.seqera.wave.plugin.cli

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.wave.plugin.WaveClient
import nextflow.Session
import nextflow.container.DockerBuilder
import nextflow.exception.AbortOperationException

/**
 * Implements Wave debug command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveRunCmd {

    private Session session

    private Map containerParams

    private List<Path> containerMounts

    private Set<String> environment = new HashSet<>()

    WaveRunCmd(Session session) { this.session=session }

    WaveRunCmd withContainerParams(Map params) {
        this.containerParams = params
        return this
    }

    WaveRunCmd withMounts(List<Path> path) {
        this.containerMounts = path
        return this
    }

    WaveRunCmd withEnvironment(String... envs) {
        this.environment.addAll(envs)
        return this
    }

    void runContainer(List<String> args) {
        if( !args )
            throw new AbortOperationException("Missing container image - usage: nextflow plugin nf-wave:run-container <image>")
        final image = args.pop()
        final target = resolveTargetImage(image)
        log.info "Resolved image: '$image' => '$target'"
        runContainer(target, args)
    }

    void runContainer(String image, List<String> args=Collections.emptyList()) {
        final containerConfig = session.getContainerConfig()
        final containerBuilder = new DockerBuilder(image)
                .addMountWorkDir(false)
                .addRunOptions('--rm')
                .addMounts(containerMounts)
                .params(containerConfig)
                .params(containerParams)

        // add env variables
        environment.addAll( containerConfig.getEnvWhitelist() )
        for( String env : environment )
            containerBuilder.addEnv(env)

        // assemble the final command
        final containerCmd = containerBuilder
                .build()
                .getRunCommand(args.join(' '))
                .replaceAll('-w "\\$NXF_TASK_WORKDIR" ','') // <-- hack to remove the PWD work dir

        log.debug "Running: $containerCmd"
        final process = new ProcessBuilder()
                .command(['sh','-c',containerCmd])
                .directory(Path.of(".").toFile())
                .inheritIO()
                .start()

        process.waitFor()
    }

    protected String resolveTargetImage(String image) {
        final resp = new WaveClient(session).sendRequest(image)
        return resp?.targetImage
    }

    void getContainer(List<String> args) {
        if( !args )
            throw new AbortOperationException("Missing container image - usage: nextflow plugin nf-wave:get-container <image>")
        final image = args.pop()
        final target = resolveTargetImage(image)
        log.info """\
                Source container: $image
                Waved  container: $target""".stripIndent(true)
    }
}
