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
 */

package io.seqera.wave.plugin.cli

import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.wave.plugin.ContainerConfig
import io.seqera.wave.plugin.packer.Packer
import io.seqera.wave.plugin.util.BasicCliOpts
import nextflow.cli.PluginCommandBase
import nextflow.exception.AbortOperationException
import nextflow.io.BucketParser

/**
 * Implements Wave CLI as a top-level command
 *
 * @author Edmund Miller <edmund@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Execute Wave container operations")
class WaveCmd extends PluginCommandBase {

    static final String COMMAND_NAME = 'wave'
    static final String COMMAND_DESCRIPTION = 'Execute Wave container operations'
    static final String PLUGIN_ID = 'nf-wave'

    @Parameter(hidden = true)
    List<String> args

    WaveCmd() {
        super(PLUGIN_ID)
    }

    @Override
    String getCommandName() {
        return COMMAND_NAME
    }

    @Override
    String getCommandDescription() {
        return COMMAND_DESCRIPTION
    }

    @Override
    protected void execute() {
        if (!args) {
            showHelp()
            return
        }

        final subCommand = args[0]
        final subArgs = args.size() > 1 ? args[1..-1] : []
        
        // Handle help requests for subcommands
        if (subCommand == 'help' || subCommand == '-h' || subCommand == '--help') {
            showHelp()
            return
        }
        
        // Execute the subcommand directly using the session
        switch (subCommand) {
            case 'get-container':
                new WaveRunCmd(getSession()).getContainer(subArgs)
                break
            case 'run-container':
                new WaveRunCmd(getSession()).runContainer(subArgs)
                break
            case 'pack':
                println packContainer(subArgs)
                break
            case 'debug-task':
                new WaveDebugCmd(getSession()).taskDebug(subArgs)
                break
            default:
                throw new AbortOperationException("Unknown wave command: $subCommand")
        }
    }

    private void showHelp() {
        println """\
Usage: nextflow wave <subcommand> [args...]

Execute Wave container operations

Available subcommands:
  get-container     Get a container image URL from Wave service
  run-container     Run a container using Wave service
  pack              Pack a directory into a Wave container layer
  debug-task        Debug a specific task execution with Wave

Examples:
  nextflow wave get-container --image ubuntu:20.04
  nextflow wave run-container --image ubuntu:20.04 --run 'echo hello'
  nextflow wave pack /path/to/directory
  nextflow wave debug-task <task-id>

For detailed help on a specific subcommand, use:
  nextflow wave <subcommand> --help"""
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

    @Override
    int getPriority() {
        return 200  // Higher than default plugin priority to take precedence
    }

}