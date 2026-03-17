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

package nextflow.cli.module

import com.beust.jcommander.IParameterValidator
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.Parameters
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleChannel
import io.seqera.npr.api.schema.v1.ModuleChannelItem
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.api.schema.v1.ModuleRelease
import io.seqera.npr.api.schema.v1.ModuleTool
import nextflow.cli.CmdBase
import nextflow.config.ConfigBuilder
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleReference
import io.seqera.npr.client.RegistryClient
import nextflow.module.RegistryClientFactory
import nextflow.util.TestOnly

import java.nio.file.Path
import java.nio.file.Paths

/**
 * Module info subcommand - displays module metadata and usage template
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Show module information and usage template")
class CmdModuleInfo extends CmdBase {

    @Parameter(names = ["-version"], description = "Module version")
    String version

    @Parameter(
        names = ['-o', '-output'],
        description = 'Output mode for reporting search results: text, json',
        validateWith = OutputModeValidator
    )
    String output = 'text'

    static class OutputModeValidator implements IParameterValidator {

        private static final List<String> MODES = List.of('text', 'json')

        @Override
        void validate(String name, String value) {
            if( !MODES.contains(value) )
                throw new ParameterException("Output mode must be one of $MODES (found: $value)")
        }
    }

    @Parameter(description = "[scope/name]", required = true)
    List<String> args

    @TestOnly
    protected Path root

    @TestOnly
    protected RegistryClient client

    @Override
    String getName() {
        return 'info'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 ) {
            throw new AbortOperationException("Incorrect number of arguments")
        }

        def reference = ModuleReference.parse(args[0])

        // Get config
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(baseDir)
            .build()
        final registryConfig = config.navigate('registry') as RegistryConfig ?: new RegistryConfig()


        // Fetch full metadata from registry to get input/output parameters
        def registryClient = this.client ?: RegistryClientFactory.forConfig(registryConfig)
        ModuleRelease release = null

        try {
            if( version ) {
                release = registryClient.fetchRelease(reference.fullName, version)

            } else {
                release = registryClient.fetchModule(reference.fullName).latest
            }
        } catch( Exception e ) {
            log.warn "Failed to fetch metadata from registry: ${e.message}"
        }
        if( !release ) {
            throw new AbortOperationException("No release information available for ${reference}")
        }
        if( !release.metadata ) {
            log.info("No metadata found for $reference ${release.version ? "($release.version)" : ''}")
        }
        def moduleUrl = buildModuleUrl(registryConfig.url, reference, release.version)
        if( !output || output == 'text' ) {
            printFormattedInfo(reference, release, moduleUrl)
        } else if( output == 'json' ) {
            printJsonInfo(reference, release, moduleUrl)
        } else {
            throw new AbortOperationException("Not implemented output mode $output)")
        }
    }

    private void printFormattedInfo(ModuleReference reference, ModuleRelease release, String moduleUrl) {
        ModuleMetadata metadata = release.metadata
        println ""
        println "Module:      ${reference}"
        println "Version:     ${release.version}"
        println "URL:         ${moduleUrl}"
        println "Description: ${metadata.description ?: release.description ?: 'N/A'}"

        if( metadata.authors ) {
            println "Authors:     ${metadata.authors.join(', ')}"
        }

        if( metadata.maintainers ) {
            println "Maintainers: ${metadata.maintainers.join(', ')}"
        }

        if( metadata.keywords ) {
            println "Keywords:    ${metadata.keywords.join(', ')}"
        }

        printToolsInfo(metadata?.tools ?: [])

        printInputsInfo(metadata.input ?: [])

        printOutputsInfo(metadata.output ?: [:])

        // Generate and display usage template
        println ""
        println "Usage Template:"
        println "-" * 80
        println generateUsageTemplate(reference, metadata).join(" \\\n    ")
        println ""
    }

    private void printOutputsInfo(Map<String, ModuleChannel> outputs) {
        if( outputs ) {
            println ""
            println "Output:"
            outputs.each { name, output ->
                println "- ${name} ${output.tuple ? '(tuple)' : ''}"
                displayChannel("\t", output)
            }
        }
    }

    private void printInputsInfo(List<ModuleChannel> inputs) {
        if( inputs ) {
            println ""
            println "Input:"
            inputs.each { input ->
                if( input.tuple ) {
                    println "- (tuple)"
                    displayChannel("\t", input)
                } else
                    displayChannel("", input)
            }
        }
    }

    private void printToolsInfo(List<ModuleTool> toolsList) {
        if( toolsList ) {
            println ""
            println "Tools:"
            toolsList.each { tool ->
                println "  - ${tool.name}${tool.version ? ' v' + tool.version : ''}"
                if( tool.homepage ) {
                    println "    Homepage: ${tool.homepage}"
                }
            }
        }
    }

    private void displayChannel(String prefix, ModuleChannel channel) {
        channel.items.each { ModuleChannelItem item ->
            println "${prefix}- ${item.name}${item.type ? ' (' + item.type + ')' : ''}"
            if( item.description ) {
                println "${prefix}\t${item.description.replaceAll(/\R/, ' ')}"
            }
            if( item.pattern ) {
                println "${prefix}\tPattern: ${item.pattern}"
            }
        }
    }

    private List<String> generateUsageTemplate(ModuleReference reference, ModuleMetadata metadata) {
        def template = new ArrayList<String>()
        template.add("nextflow module run ${reference}".toString())
        if( version )
            template.add(" -version $version".toString())

        def inputs = metadata?.input ?: []
        inputs.each { input ->
            input.items.each { ModuleChannelItem item ->
                template.add(reference.scope == 'nf-core'
                    ? inferNfCoreParam(item.name, item.type)
                    : inferNormalParam(item.name, item.type))

            }
        }
        if( reference.scope == 'nf-core' ) {
            template.add('--outdir <OUTPUT_DIRECTORY>')
        }
        return template
    }

    private static String inferNfCoreParam(String paramName, String type) {
        if( type?.equalsIgnoreCase("map") && paramName.equalsIgnoreCase("meta") ) {
            return "--${paramName}.id <ID>"
        }
        if( type?.equalsIgnoreCase("file") || type?.equalsIgnoreCase("path") ) {
            return "--${paramName} ${inferBioFilePlaceholder(paramName)}"
        }
        return inferNormalParam(paramName, type)
    }

    private static String inferBioFilePlaceholder(String paramName) {
        final String lower = paramName.toLowerCase()
        if( lower.contains("fasta") ) return "<FASTA_FILE>"
        if( lower.contains("bam") ) return "<BAM_FILE>"
        if( lower.contains("fastq") || lower.equals("reads") ) return "<FASTQ_FILE>"
        if( lower.contains("vcf") ) return "<VCF_FILE>"
        if( lower.contains("ref") ) return "<REFERENCE_FILE>"
        if( lower.contains("bed") ) return "<BED_FILE>"
        if( lower.contains("gff") || lower.contains("gtf") ) return "<ANNOTATION_FILE>"

        return "<${paramName.toUpperCase().replaceAll(/[^A-Z0-9]/, '_')}_PATH>"
    }

    private static String inferNormalParam(String paramName, String type) {
        final paramPlaceholder = paramName.toUpperCase().replaceAll(/[^A-Z0-9]/, '_')
        if( type?.equalsIgnoreCase("map") ) {
            return "--${paramName}.<KEY> <${paramPlaceholder}_KEY_VALUE>"
        }
        if( type?.equalsIgnoreCase("file") || type?.equalsIgnoreCase("path") ) {
            return "--${paramName} <${paramPlaceholder}_PATH>"
        }
        return "--${paramName} <${paramPlaceholder}>"
    }

    private static String buildModuleUrl(String registryUrl, ModuleReference reference, String version) {
        // Strip /api suffix to get the base UI URL
        def baseUrl = registryUrl.endsWith('/api') ? registryUrl[0..-5] : registryUrl
        def encodedName = URLEncoder.encode(reference.name, 'UTF-8')
        return "${baseUrl}/admin/modules/${reference.scope}/${encodedName}@${version}"
    }

    private void printJsonInfo(ModuleReference reference, ModuleRelease release, String moduleUrl) {
        def metadata = release?.metadata
        def info = [
            name       : reference.toString(),
            fullName   : reference.fullName,
            version    : release.version,
            url        : moduleUrl,
            description: metadata.description ?: release.description,
            authors    : metadata.authors,
            keywords   : metadata.keywords,
        ]

        def toolsList = metadata.tools ?: []
        if( toolsList ) {
            info.tools = toolsList.collect { tool ->
                return [
                    name         : tool.name,
                    version      : tool.version,
                    homepage     : tool.homepage,
                    documentation: tool.documentation
                ]
            }
        }

        def inputs = metadata.input ?: []
        if( inputs ) {
            info.input = inputs.collect { input ->
                return [
                    tuple: input.tuple,
                    items: input.items?.collect { item ->
                        [
                            name       : item.name,
                            type       : item.type,
                            description: item.description,
                            pattern    : item.pattern
                        ]
                    }
                ]
            }
        }

        def outputs = metadata.output ?: [:]
        if( outputs ) {
            info.output = outputs.collectEntries { name, output ->
                return [name, [
                    tuple: output.tuple,
                    items: output.items?.collect { item ->
                        [
                            name       : item.name,
                            type       : item.type,
                            description: item.description,
                            pattern    : item.pattern
                        ]
                    }
                ]]
            }
        }

        info.usageTemplate = generateUsageTemplate(reference, metadata).join(" ")

        println JsonOutput.prettyPrint(JsonOutput.toJson(info))
    }
}
