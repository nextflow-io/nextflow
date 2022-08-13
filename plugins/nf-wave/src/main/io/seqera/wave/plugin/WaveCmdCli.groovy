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

package io.seqera.wave.plugin

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.PluginAbstractExec
import nextflow.container.DockerBuilder
import nextflow.exception.AbortOperationException
/**
 * Implements Wave CLI tool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveCmdCli implements PluginAbstractExec {

    List<String> getCommands() {  ['get-container','run-container'] }
    @Override
    int exec(String cmd, List<String> args) {

        if( cmd == 'get-container' )
            getContainer(args)
        else
            runContainer(args)

        return 0
    }

    protected void getContainer(List<String> args) {
        if( !args )
            throw new AbortOperationException("Missing container image - usage: nextflow plugin exec nf-wave get-container <image>")
        final image = args.pop()
        final target = resolveTargetImage(image)
        log.info """\
                Source container: $image
                Waved  container: $target""".stripIndent()
    }

    protected void runContainer(List<String> args) {
        if( !args )
            throw new AbortOperationException("Missing container image - usage: nextflow plugin exec nf-wave container-run <image>")
        final image = args.pop()
        final target = resolveTargetImage(image)
        log.info "Resolved image: '$image' => '$target'"
        final containerConfig = getSession().getContainerConfig()
        final containerBuilder = new DockerBuilder(target)
                .params(containerConfig)
        // add env variables
        for( String env : containerConfig.getEnvWhitelist())
            containerBuilder.addEnv(env)
        // assemble the final command
        final containerCmd = containerBuilder
                .build()
                .getRunCommand(args.join(' '))

        log.info "Running: $containerCmd"
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


}
