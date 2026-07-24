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

import java.nio.file.Files
import java.nio.file.Path

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleInfo

/**
 * Module create subcommand -- creates a new module skeleton
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Create a new module skeleton")
class CmdModuleCreate extends CmdBase {

    @Parameter(description = "[namespace/name]")
    List<String> args

    @Parameter(names = ['-kind'], description = "Module kind: Process (default) or Workflow")
    String kind

    @Parameter(names = ['-typed'], description = "Generate a statically-typed module (usable directly with `nextflow module run`)", arity = 0)
    boolean typed

    @Override
    String getName() {
        return 'create'
    }

    @Override
    void run() {
        String namespace
        String name

        if( args && args.size() == 1 && args[0].contains('/') ) {
            // non-interactive: namespace/name passed as argument
            final slash = args[0].indexOf('/')
            namespace = args[0].substring(0, slash)
            name = args[0].substring(slash + 1)
            if( !namespace || !name )
                throw new AbortOperationException("Invalid module identifier -- expected format: namespace/name")
        }
        else if( args ) {
            throw new AbortOperationException("Invalid arguments -- usage: nextflow module create [namespace/name]")
        }
        else {
            // interactive mode
            print "Enter module namespace: "
            namespace = readLine()?.trim()
            if( !namespace )
                throw new AbortOperationException("Module namespace cannot be empty")

            print "Enter module name: "
            name = readLine()?.trim()
            if( !name )
                throw new AbortOperationException("Module name cannot be empty")

            println ""
            println "  Module namespace : $namespace"
            println "  Module name      : $name"
            println "  Directory        : ./modules/$namespace/$name"
            println ""
            print "Are you OK to continue [y/N]? "
            final confirm = readLine()
            if( confirm?.toLowerCase() != 'y' ) {
                println "Module creation aborted."
                return
            }
        }

        validateSegment('namespace', namespace)
        validateSegments('name', name)
        createModule(namespace, name, normalizeKind(kind), typed)
    }

    static private String normalizeKind(String value) {
        if( !value )
            return 'Process'
        final k = value.toLowerCase().capitalize()
        if( k != 'Process' && k != 'Workflow' )
            throw new AbortOperationException("Invalid module kind '${value}' -- must be 'Process' or 'Workflow'")
        return k
    }

    static private void validateSegment(String field, String value) {
        if( !value.matches('[a-zA-Z0-9][a-zA-Z0-9._\\-]*') )
                throw new AbortOperationException("Invalid module $field '${value}' -- only alphanumeric characters, hyphens, underscores and dots are allowed, and must start with an alphanumeric character")
    }

    static private void validateSegments(String field, String value) {
        for( String segment : value.tokenize('/') ) {
            validateSegment(field, segment)
        }
    }

    protected Path modulesBase() {
        return Path.of('modules')
    }

    protected void createModule(String namespace, String name, String kind = 'Process', boolean typed = false) {
        final moduleDir = modulesBase().resolve(namespace).resolve(name)
        if( Files.exists(moduleDir) )
            throw new AbortOperationException("Module directory already exists: $moduleDir")

        // create directory structure
        Files.createDirectories(moduleDir)

        // create main.nf
        moduleDir.resolve('main.nf').text = mainNf(namespace, name, kind, typed)

        // create README.md
        moduleDir.resolve('README.md').text = readmeMd(namespace, name)

        // create meta.yml
        moduleDir.resolve('meta.yml').text = metaYml(namespace, name, kind, typed)

        // create .module-info so it's recognised as a Nextflow managed module
        Files.createFile(moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE))

        final defName = name.replaceAll('[^a-zA-Z0-9_]', '_').toUpperCase()
        println "Module created successfully at path: $moduleDir"
        println ""
        // an untyped workflow module cannot be run directly (`module run` requires typed take:/emit:)
        if( kind == 'Workflow' && !typed ) {
            println "Include the workflow module in a pipeline:"
            println ""
            println "  include { $defName } from '$namespace/$name'"
        }
        else {
            println "To run the module:"
            println ""
            println "  nextflow module run $namespace/$name --greeting 'Hello world!'"
        }
    }

    static private String readLine() {
        final console = System.console()
        return console != null
            ? console.readLine()
            : new BufferedReader(new InputStreamReader(System.in)).readLine()
    }

    static String mainNf(String namespace, String name, String kind = 'Process', boolean typed = false) {
        final defName = name.replaceAll('[^a-zA-Z0-9_]', '_').toUpperCase()

        if( kind == 'Workflow' && typed ) {
            return """\
            /*
             * Workflow module: ${namespace}/${name}
             * TODO: rename the workflow, replace the example take/emit and types, and implement the logic.
             */

            nextflow.enable.types = true

            workflow ${defName} {
                take:
                greeting: String

                main:
                // TODO: implement the workflow logic
                message = greeting

                emit:
                result: String = message
            }
            """.stripIndent()
        }

        if( kind == 'Workflow' ) {
            return """\
            /*
             * Workflow module: ${namespace}/${name}
             * TODO: rename the workflow, replace the example take/emit, and implement the logic.
             */

            workflow ${defName} {
                take:
                ch_input

                main:
                // TODO: implement the workflow logic
                ch_output = ch_input

                emit:
                output = ch_output
            }
            """.stripIndent()
        }

        if( typed ) {
            return """\
            /*
             * Module: ${namespace}/${name}
             * TODO: rename the process, replace the example input/output and types, and implement the script.
             */

            nextflow.enable.types = true

            process ${defName} {
                input:
                greeting: String

                output:
                message: String = stdout()

                script:
                \"\"\"
                echo '\${greeting}'
                \"\"\"
            }
            """.stripIndent()
        }

        return """\
        /*
         * Module: ${namespace}/${name}
         * TODO: rename the process, replace the example input/output, and implement the script.
         */

        process ${defName} {
            input:
            val greeting

            output:
            stdout

            script:
            \"\"\"
            echo '\${greeting}'
            \"\"\"
        }
        """.stripIndent()
    }

    static String metaYml(String namespace, String name, String kind = 'Process', boolean typed = false) {
        if( kind == 'Workflow' && typed ) {
            // typed workflow: derive input/output from the scaffold's take:/emit:
            // (typed workflows require Nextflow 26.04.0+)
            return """\
            name: ${namespace}/${name}
            version: 1.0.0
            kind: Workflow
            description: A brief description of the ${namespace}/${name} workflow module
            license: Apache-2.0
            requires:
              nextflow: ">=26.04.0"
            input:
              - name: greeting
                type: string
                description: A greeting string
            output:
              - name: result
                type: string
                description: The greeting message
            """.stripIndent()
        }
        if( kind == 'Workflow' ) {
            // untyped workflow: take/emit have no declared types, but the generated scaffold's
            // take/emit are channels, so the interface is documented with the generic channel type
            return """\
            name: ${namespace}/${name}
            version: 1.0.0
            kind: Workflow
            description: A brief description of the ${namespace}/${name} workflow module
            license: Apache-2.0
            requires:
              nextflow: ">=24.04.0"
            input:
              - name: ch_input
                type: channel
                description: The input channel
            output:
              - name: output
                type: channel
                description: The output channel
            """.stripIndent()
        }
        if( typed ) {
            // typed process (requires Nextflow 25.10.0+)
            return """\
            name: ${namespace}/${name}
            version: 1.0.0
            description: A brief description of the ${namespace}/${name} module
            license: Apache-2.0
            requires:
              nextflow: ">=25.10.0"
            input:
              - name: greeting
                type: string
                description: A greeting string
            output:
              - name: message
                type: string
                description: The greeting message
            """.stripIndent()
        }
        return """\
        name: ${namespace}/${name}
        version: 1.0.0
        description: A brief description of the ${namespace}/${name} module
        license: Apache-2.0
        input:
          - name: greeting
            type: string
            description: A greeting string
        output:
          - name: stdout
            type: string
            description: The greeting message
        """.stripIndent()
    }

    static String readmeMd(String namespace, String name) {
        """\
        # ${namespace}/${name}

        ## Summary

        A brief description of the `${namespace}/${name}` module.

        ## Get started

        Include this module in your Nextflow pipeline:

        ```nextflow
        include { ${name.replaceAll('[^a-zA-Z0-9_]', '_').toUpperCase()} } from '${namespace}/${name}'
        ```

        ## Dependencies

        None.

        ## License

        Apache-2.0
        """.stripIndent()
    }
}
