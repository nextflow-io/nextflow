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
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleSearchResult
import io.seqera.npr.api.schema.v1.PipelineSummary
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.config.ConfigCmdAdapter
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleSpec
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
@Parameters(commandDescription = "Search for modules and pipelines in the registry")
class CmdModuleSearch extends CmdBase {

    private static final int MAX_PAGE_SIZE = 100

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
        def config = new ConfigCmdAdapter()
            .setOptions(launcher.options)
            .setBaseDir(baseDir)
            .build()

        final registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()

        // Create client to search
        final client = this.client ?: RegistryClientFactory.forConfig(registryConfig)

        try {
            println "Searching for '${query}'..."
            // modules first, pipelines fill the rest of the limit
            final results = client.searchModules(query, limit)
            final modules = (results?.results ?: List.<ModuleSearchResult>of()).take(limit)
            final pipelines = searchPipelines(client, query, limit - modules.size())

            if( !modules && !pipelines.pipelines ) {
                println "No modules or pipelines found"
                return
            }

            if( !output || output == 'simple' ) {
                printFormattedResults(modules, pipelines.pipelines)
            } else if( output == 'json' ) {
                // the registry totals; pipelines are counted only when listed
                final total = (results?.totalResults ?: modules.size()) + (pipelines.total ?: pipelines.pipelines.size())
                printJsonResults(query, total, modules, pipelines.pipelines)
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

    /**
     * The module kind, defaulting to {@code Process} when the registry does not
     * report one (older registries or process modules that predate the field).
     */
    private static String kindOf(ModuleSearchResult result) {
        return result.kind != null ? result.kind.toString() : 'Process'
    }

    @Canonical
    private static class PipelineResults {
        List<PipelineSummary> pipelines = []
        Long total
    }

    /**
     * At most {@code max} pipelines matching the query, with the registry total;
     * empty when the registry can't list pipelines
     */
    private static PipelineResults searchPipelines(RegistryClient client, String query, int max) {
        final result = new PipelineResults()
        if( max <= 0 )
            return result
        // the registry caps the page size, and the size must be the same on every page
        final size = Math.min(max, MAX_PAGE_SIZE)
        try {
            for( int page = 1; result.pipelines.size() < max; page++ ) {
                final response = client.listPipelines(query, page, size)
                if( result.total == null )
                    result.total = response?.totalSize
                if( !response?.pipelines )
                    break
                result.pipelines.addAll(response.pipelines)
                if( response.last )
                    break
            }
        }
        catch( Exception e ) {
            log.debug "Failed to search pipelines: ${e.message}"
        }
        result.pipelines = result.pipelines.take(max)
        return result
    }

    private void printFormattedResults(List<ModuleSearchResult> modules, List<PipelineSummary> pipelines) {
        println ""
        println "Top ${modules.size() + pipelines.size()} matching result(s):"
        println ""

        modules.each { ModuleSearchResult result ->
            println "  ${result.name}"
            println "    Kind: ${kindOf(result)}"
            if( result.description ) {
                println "    Description: ${result.description}"
            }
            println ""
        }
        pipelines.each { PipelineSummary p ->
            println "  ${p.name}"
            println "    Kind: ${ModuleSpec.KIND_PIPELINE}"
            println "    Version: ${p.latestVersion}"
            if( p.metadata?.description ) {
                println "    Description: ${p.metadata.description}"
            }
            println ""
        }
    }

    private void printJsonResults(String query, long total, List<ModuleSearchResult> modules, List<PipelineSummary> pipelines) {
        final rows = modules.collect { ModuleSearchResult result ->
            [
                name          : result.name,
                repositoryPath: result.repositoryPath,
                kind          : kindOf(result),
                description   : result.description,
                relevanceScore: result.relevanceScore,
                keywords      : result.keywords,
                tools         : result.tools,
                revoked       : result.revoked
            ]
        } as List<Map>
        pipelines.each { PipelineSummary p ->
            rows << [
                name         : p.name,
                kind         : ModuleSpec.KIND_PIPELINE,
                latestVersion: p.latestVersion,
                description  : p.metadata?.description,
                keywords     : p.metadata?.keywords
            ]
        }

        println JsonOutput.prettyPrint(JsonOutput.toJson(
            query: query,
            totalResults: total,
            results: rows
        ))
    }
}
