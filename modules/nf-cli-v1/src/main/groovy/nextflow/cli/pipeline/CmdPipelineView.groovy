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

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import io.seqera.npr.client.RegistryClient
import nextflow.cli.CmdBase
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.RegistryClientFactory
import nextflow.pipeline.PipelineManager
import nextflow.util.TestOnly

import java.nio.file.Paths

/**
 * Pipeline view subcommand - displays pipeline info and parameters
 */
@CompileStatic
@Parameters(commandDescription = "Show pipeline information")
class CmdPipelineView extends CmdBase {

    @Parameter(names = ["-version"], description = "Pipeline version (default: latest)")
    String version

    @Parameter(
        names = ['-o', '-output'],
        description = PipelineCmdSupport.OUTPUT_MODE_HELP,
        validateWith = PipelineCmdSupport.OutputModeValidator
    )
    String output = 'text'

    @Parameter(description = "<scope/name>", required = true)
    List<String> args

    @TestOnly
    protected RegistryClient client

    @Override
    String getName() {
        return 'view'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 )
            throw new AbortOperationException("Unexpected number of parameters")

        final name = args[0]
        final registryConfig = PipelineCmdSupport.registryConfig(this, Paths.get('.').toAbsolutePath().normalize())
        final client = this.client ?: RegistryClientFactory.forConfig(registryConfig)
        final pipeline = client.getPipeline(name)
        final release = version ? client.getPipelineRelease(name, version) : pipeline?.latest
        if( !pipeline || !release )
            throw new AbortOperationException("No release information available for ${name}")

        final url = pipelineUrl(registryConfig.url, name, release.version)
        final params = schemaParams(release.schema)
        final usage = "nextflow run ${name} -hub ${URI.create(registryConfig.url).host} -r ${release.version}".toString()

        if( output == 'json' ) {
            final info = [
                name         : name,
                version      : release.version,
                repositoryUrl: pipeline.repositoryUrl,
                gitRevision  : release.gitRevision,
                commitId     : release.commitId,
                url          : url,
                description  : pipeline.metadata?.description,
                keywords     : pipeline.metadata?.keywords,
                params       : params,
                usage        : usage
            ]
            println JsonOutput.prettyPrint(JsonOutput.toJson(info))
            return
        }

        println "Name:        ${name}"
        println "Version:     ${release.version}"
        println "Repository:  ${pipeline.repositoryUrl}"
        println "Revision:    ${release.gitRevision}"
        println "Commit:      ${release.commitId}"
        println "URL:         ${url}"
        if( pipeline.metadata?.description )
            println "Description: ${pipeline.metadata.description}"
        if( pipeline.metadata?.keywords )
            println "Keywords:    ${pipeline.metadata.keywords.join(', ')}"
        if( params ) {
            println ""
            println "Parameters:"
            params.each { p ->
                println "  --${p.name}${p.type ? ' (' + p.type + ')' : ''}"
                if( p.description )
                    println "      ${p.description}"
            }
        }
        println ""
        println "Usage:"
        println "  ${usage}"
    }

    /** Lists params from `$defs`/`definitions` groups and top-level `properties`; types may be missing. */
    protected static List<Map<String,String>> schemaParams(Map schema) {
        final result = new ArrayList<Map<String,String>>()
        if( !schema )
            return result
        final groups = (schema['$defs'] ?: schema['definitions']) as Map<String,Object>
        final props = new ArrayList<Map>()
        groups?.values()?.each { g ->
            final p = (g as Map)?.get('properties')
            if( p instanceof Map )
                props << (p as Map)
        }
        if( schema.get('properties') instanceof Map )
            props << (schema.get('properties') as Map)
        props.each { Map m ->
            m.each { k, v ->
                final d = v instanceof Map ? (Map) v : [:]
                result << [
                    name       : k as String,
                    type       : d['type'] as String,
                    description: (d['description'] as String)?.replaceAll(/\R/, ' ')
                ]
            }
        }
        return result
    }

    private static String pipelineUrl(String registryUrl, String name, String version) {
        return "${PipelineManager.pageUrl(registryUrl)}/pipelines/${name}@${version}"
    }
}
