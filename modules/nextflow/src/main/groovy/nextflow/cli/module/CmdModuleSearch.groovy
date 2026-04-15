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

package nextflow.cli.module

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleSearchResult
import io.seqera.npr.api.schema.v1.SearchModulesResponse
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.util.TestOnly

import java.nio.file.Paths

/**
 * Module search subcommand
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Search for modules in the registry")
class CmdModuleSearch extends CmdBase {

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
        if( !args || args.size() != 1 ) {
            throw new AbortOperationException("Unexpected number of parameters")
        }
        String query = args[0]

        // Get config
        def baseDir = Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(baseDir)
            .build()

        final registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()

        // Create client to search
        final client = this.client ?: RegistryClientFactory.forConfig(registryConfig)

        try {
            println "Searching for '${query}'..."
            final results = client.searchModules(query, limit)

            if( !results || results.totalResults == 0 || !results.results || results.results.isEmpty() ) {
                println "No modules found"
                return
            }

            if( !output || output == 'simple' ) {
                printFormattedResults(results)
            } else if( output == 'json' ) {
                printJsonResults(results)
            } else {
                throw new AbortOperationException("Not implemented output mode $output)")
            }
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            log.error("Failed to search modules", e)
            throw new AbortOperationException("Search failed: ${e.message}", e)
        }
    }

    private void printFormattedResults(SearchModulesResponse response) {
        println ""
        println "Top ${response.totalResults} matching module(s):"
        println ""

        response.results.each { ModuleSearchResult result ->
            println "  ${result.name}"
            if( result.description ) {
                println "    Description: ${result.description}"
            }
            println ""
        }
    }

    private void printJsonResults(SearchModulesResponse response) {
        final modules = response.results.collect { ModuleSearchResult result ->
            [
                name          : result.name,
                repositoryPath: result.repositoryPath,
                description   : result.description,
                relevanceScore: result.relevanceScore,
                keywords      : result.keywords,
                tools         : result.tools,
                revoked       : result.revoked
            ]
        }

        println JsonOutput.prettyPrint(JsonOutput.toJson(
            query: response.query,
            totalResults: response.totalResults,
            results: modules
        ))
    }
}
