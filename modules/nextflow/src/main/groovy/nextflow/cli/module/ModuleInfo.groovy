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

import com.beust.jcommander.Parameter
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
import nextflow.module.InstalledModule
import nextflow.module.ModuleReference
import nextflow.module.ModuleRegistryClient
import nextflow.module.ModuleSpec
import nextflow.module.ModuleStorage
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
class ModuleInfo extends CmdBase {

    @Parameter(names = ["-version"], description = "Module version")
    String version

    @Parameter(names = ["-json"], description = "Output in JSON format", arity = 0)
    boolean jsonOutput = false

    @Parameter(description = "[scope/name]", required = true)
    List<String> args

    @TestOnly
    protected Path root

    @TestOnly
    protected ModuleRegistryClient client

    @Override
    String getName() {
        return 'info'
    }

    @Override
    void run() {
        if( !args || args.size() != 1 ) {
            throw new AbortOperationException("Incorrect number of arguments")
        }

        def moduleRef = '@' + args[0]
        def reference = ModuleReference.parse(moduleRef)

        // Get config
        def baseDir = root ?: Paths.get('.').toAbsolutePath().normalize()
        def config = new ConfigBuilder()
            .setOptions(launcher.options)
            .setBaseDir(baseDir)
            .build()
        final registryConfig = config.navigate('registry') as RegistryConfig


        // Fetch full metadata from registry to get input/output parameters
        def registryClient = this.client ?: new ModuleRegistryClient(registryConfig)
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
        if( release?.metadata ) {
            log.info("No metadata found for $reference.nameWithoutPrefix ${release?.version ? "($release.version)" : ''}")
        }
        if( jsonOutput ) {
            printJsonInfo(reference, release)
        } else {
            printFormattedInfo(reference, release)
        }
    }

    private void printFormattedInfo(ModuleReference reference, ModuleRelease release) {
        ModuleMetadata metadata = release.metadata
        println ""
        println "Module:      ${reference.nameWithoutPrefix}"
        println "Version:     ${release.version}"
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
        println generateUsageTemplate(reference, metadata)
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

    private String generateUsageTemplate(ModuleReference reference, ModuleMetadata metadata) {
        def template = new StringBuilder()
        template.append("nextflow module run ${reference.nameWithoutPrefix}")
        if( version )
            template.append(" -version $version")

        // Use metadata from registry if available, otherwise use spec from meta.yml
        def inputs = metadata?.input ?: []

        if( inputs ) {
            inputs.each { input ->
                input.items.each { ModuleChannelItem item ->
                    def placeholder = item.name.toUpperCase().replaceAll(/[^A-Z0-9_]/, '_')
                    if( item.type.equalsIgnoreCase("map") ) {
                        template.append(" --${item.name}.<key> <${placeholder}_KEY>")
                    } else {
                        template.append(" --${item.name} <${placeholder}>")
                    }
                }
            }
        }
        return template.toString()
    }

    private void printJsonInfo(ModuleReference reference, ModuleRelease release) {
        def metadata = release?.metadata
        def info = [
            name       : reference.nameWithoutPrefix,
            fullName   : reference.fullName,
            version    : release.version,
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

        info.usageTemplate = generateUsageTemplate(reference, metadata)

        println JsonOutput.prettyPrint(JsonOutput.toJson(info))
    }
}