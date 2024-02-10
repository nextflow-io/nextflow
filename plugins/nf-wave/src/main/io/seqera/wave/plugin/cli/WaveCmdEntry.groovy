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

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.wave.plugin.ContainerConfig
import io.seqera.wave.plugin.packer.Packer
import io.seqera.wave.plugin.util.BasicCliOpts
import nextflow.cli.PluginAbstractExec
import nextflow.exception.AbortOperationException
import nextflow.io.BucketParser
/**
 * Implements Wave CLI tool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WaveCmdEntry implements PluginAbstractExec {

    List<String> getCommands() {  ['get-container','run-container', 'debug-task', 'pack'] }

    @Override
    int exec(String cmd, List<String> args) {

        switch (cmd) {
            case 'get-container':
                new WaveRunCmd(session).getContainer(args)
                break
            case 'run-container':
                new WaveRunCmd(session).runContainer(args)
                break
            case 'pack':
                println packContainer(args)
                break
            case 'debug-task':
                new WaveDebugCmd(session).taskDebug(args)
                break
            default:
                throw new AbortOperationException("Unknown wave command: $cmd")
        }

        return 0
    }

    protected String packContainer(List<String> args) {
        final cli = BasicCliOpts.parse(args)
        final packer = new Packer()
        if( !cli.args )
            throw new AbortOperationException("Missing pack target directory")
        final root = Path.of(cli.args.pop())
        if( !Files.exists(root) )
            throw new AbortOperationException("Pack target path does not exist: $root")
        if( !Files.isDirectory(root) )
            throw new AbortOperationException("Pack target path is not a directory: $root")
        // determine target location form CLI option
        final location = cli.options.location
        // pack the layer
        final containerConfig = new ContainerConfig()
        final layer = packer.createContainerPack(root, baseName(location))
        containerConfig.appendLayer(layer)
        // set the target location
        if( location )
            layer.location = location
        // set entrypoint
        if( cli.options.entrypoint )
            containerConfig.entrypoint = [cli.options.entrypoint]
        // set work dir
        if( cli.options.workingDir )
            containerConfig.workingDir = cli.options.workingDir
        // set entrypoint
        if( cli.options.cmd )
            containerConfig.cmd = [cli.options.cmd]
        // render final object
        return JsonOutput.prettyPrint(JsonOutput.toJson(containerConfig))
    }

    static protected String baseName(String path) {
        if( !path )
            return null
        def name = BucketParser.from(path).getPath().getName()
        if( !name )
            return null
        return name
                .replaceAll(/\.gz$/,'')
                .replaceAll(/\.gzip$/,'')
                .replaceAll(/\.tar$/,'')
                .replaceAll(/\.json$/,'')
    }

}
