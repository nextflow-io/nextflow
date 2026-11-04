/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cli.pipeline

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.util.TestOnly

import java.nio.file.Paths

/**
 * Pipeline search subcommand
 */
@CompileStatic
@Parameters(commandDescription = "Search for pipelines in the registry")
class CmdPipelineSearch extends CmdBase {

    @Parameter(names = ["-limit"], description = "Maximum number of results")
    int limit = 20

    @Parameter(
        names = ['-o', '-output'],
        description = 'Output mode for reporting search results: simple, json',
        validateWith = OutputModeValidator
    )
    String output = 'simple'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('simple', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @Parameter(description = "<query>", required = true)
    List<String> args

    @TestOnly
    protected RegistryClient client

    @Override
    String getName() {
        return 'search'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")

        final client = this.client ?: RegistryClientFactory.forConfig(registryConfig())
        if( output != 'json' )
            println "Searching for '${args[0]}'..."
        final pipelines = client.searchPipelines(args[0], limit)?.pipelines ?: []

        if( output == 'json' ) {
            final rows = pipelines.collect { p ->
                [
                    name         : p.name,
                    latestVersion: p.latestVersion,
                    description  : p.metadata?.description,
                    keywords     : p.metadata?.keywords
                ]
            }
            println JsonOutput.prettyPrint(JsonOutput.toJson(rows))
            return
        }
        if( !pipelines ) {
            println "No pipelines found"
            return
        }
        println ""
        println "Top ${pipelines.size()} matching pipeline(s):"
        println ""
        pipelines.each { p ->
            println "  ${p.name}"
            println "    Version: ${p.latestVersion}"
            if( p.metadata?.description )
                println "    Description: ${p.metadata.description}"
            println ""
        }
    }

    private RegistryConfig registryConfig() {
        final config = new ConfigCmdAdapter()
            .setOptions(launcher.options)
            .setBaseDir(Paths.get('.').toAbsolutePath().normalize())
            .build()
        return config.navigate('registry') as RegistryConfig ?: new RegistryConfig()
    }
}
