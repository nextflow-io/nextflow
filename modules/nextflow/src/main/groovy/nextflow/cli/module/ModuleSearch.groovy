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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleSearchResult
import io.seqera.npr.api.schema.v1.SearchModulesResponse
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleRegistryClient
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
class ModuleSearch extends CmdBase {

    @Parameter(names = ["-limit"], description = "Maximum number of results")
    int limit = 20

    @Parameter(names = ["-json"], description = "Output in JSON format", arity=0)
    boolean jsonOutput = false

    @Parameter(description = "<query>", required = true)
    List<String> args

    @TestOnly
    protected ModuleRegistryClient client

    @Override
    String getName() {
        return 'search'
    }

    @Override
    void run() {
        if (!args && args.size() != 1 ) {
            throw new AbortOperationException("Unexpected number of parameters")
        }
        String query = args[0]

        // Get config
        def baseDir = Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(baseDir)
                .build()

        final registryConfig = config.navigate('registry') as RegistryConfig

        // Create client to search
        final client = this.client ?: new ModuleRegistryClient(registryConfig)

        try {
            println "Searching for '${query}'..."
            final results = client.search(query, limit)

            if (!results || results.totalResults == 0 || !results.results || results.results.isEmpty()) {
                println "No modules found"
                return
            }

            if (jsonOutput) {
                printJsonResults(results)
            } else {
                printFormattedResults(results)
            }
        }
        catch (AbortOperationException e) {
            throw e
        }
        catch (Exception e) {
            log.error("Failed to search modules", e)
            throw new AbortOperationException("Search failed: ${e.message}", e)
        }
    }

    private void printFormattedResults(SearchModulesResponse response) {
        println ""
        println "Found ${response.totalResults} module(s):"
        println ""

        response.results.each { ModuleSearchResult result ->
            println "  ${result.name}"
            if (result.relevanceScore != null) {
                println "    Relevance:   ${String.format('%.2f', result.relevanceScore)}"
            }
            if (result.description) {
                println "    Description: ${result.description}"
            }
            if (result.keywords && !result.keywords.isEmpty()) {
                println "    Keywords:    ${result.keywords.join(', ')}"
            }
            if (result.tools && !result.tools.isEmpty()) {
                println "    Tools:       ${result.tools.join(', ')}"
            }
            println ""
        }
    }

    private void printJsonResults(SearchModulesResponse response) {
        final modules = response.results.collect { ModuleSearchResult result ->
            [
                name: result.name,
                repositoryPath: result.repositoryPath,
                description: result.description,
                relevanceScore: result.relevanceScore,
                keywords: result.keywords,
                tools: result.tools,
                revoked: result.revoked
            ]
        }

        println JsonOutput.toJson(
            query: response.query,
            totalResults: response.totalResults,
            results: modules
        )
    }
}
