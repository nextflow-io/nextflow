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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.script.ast.ProcessNodeV1
import nextflow.script.ast.ProcessNodeV2
import nextflow.script.ast.ScriptNode
import nextflow.script.control.ScriptParser
import org.yaml.snakeyaml.Yaml

import static nextflow.module.ModuleSpec.ModuleParam

/**
 * Factory methods for module specs.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleSpecFactory {

    private static final List<String> MODULE_SPEC_FIELDS = [
        'name',
        'version',
        'description',
        'keywords',
        'license',
        'authors',
        'maintainers',
        'requires',
        'tools',
        'input',
        'output',
        'topics'
    ]

    private static final List<String> PARAM_SPEC_FIELDS = [
        'name',
        'type',
        'description',
    ]

    /**
     * Generate a module spec from the module script.
     *
     * The module name is inferred from the process name. The
     * module input/output spec is inferred from the process
     * input/output declarations.
     *
     * Additional overrides can be provided as named arguments:
     *
     * - namespace: String
     * - name: String
     * - version: String
     * - description: String
     * - license: String
     * - authors: List<String>
     *
     * @param opts
     * @param path
     * @param oldSpec
     */
    static ModuleSpec fromScript(Map opts = [:], Path path, ModuleSpec oldSpec) {
        // create initial module spec
        final spec = new ModuleSpec()
        spec.version = opts.version as String ?: oldSpec.version
        spec.description = opts.description as String ?: oldSpec.description
        spec.keywords = oldSpec.keywords
        spec.license = opts.license as String ?: oldSpec.license
        spec.authors = opts.authors as List<String> ?: oldSpec.authors
        spec.maintainers = oldSpec.maintainers
        spec.requires = oldSpec.requires
        spec.tools = oldSpec.tools
        spec._passthrough = oldSpec._passthrough

        // load script
        final parser = new ScriptParser()
        final sourceUnit = parser.parse(path.toFile())
        parser.analyze()

        final scriptNode = sourceUnit.getAST()
        if( scriptNode !instanceof ScriptNode )
            throw new AbortOperationException("Error parsing module script -- run `nextflow lint ${path}` to check for errors")

        final errors = sourceUnit.getErrorCollector().getErrors()
        if( errors != null && !errors.isEmpty() )
            throw new AbortOperationException("Error parsing module script -- run `nextflow lint ${path}` to check for errors")

        // get process definition in script
        final processes = ((ScriptNode) scriptNode).getProcesses()
        if( processes.isEmpty() )
            throw new AbortOperationException("Module script does not define any processes: ${path}")
        if( processes.size() > 1 )
            throw new AbortOperationException("Module script defines multiple processes: ${path}")

        // infer module spec properties from process definition
        final process = processes[0]

        spec.name = opts.name ?: "${opts.namespace}/${process.name.toLowerCase()}"

        if( process instanceof ProcessNodeV1 ) {
            final visitor = new ModuleSpecVisitorV1(oldSpec)
            spec.inputs = visitor.visitInputs(process)
            spec.outputs = visitor.visitOutputs(process)
            spec.topics = visitor.visitTopics(process)
        }
        else if( process instanceof ProcessNodeV2 ) {
            final visitor = new ModuleSpecVisitorV2(oldSpec)
            spec.inputs = visitor.visitInputs(process)
            spec.outputs = visitor.visitOutputs(process)
            spec.topics = visitor.visitTopics(process)
        }

        return spec
    }

    static ModuleSpec fromScript(Map opts = [:], Path path) {
        return fromScript(opts, path, new ModuleSpec())
    }

    /**
     * Load a module spec from a yaml file
     *
     * @param path
     */
    static ModuleSpec fromYaml(Path path) {
        if( !Files.exists(path) ) {
            throw new AbortOperationException("Module spec not found: ${path}")
        }

        try( final stream = Files.newInputStream(path) ) {
            final data = new Yaml().load(stream) as Map<String, Object>

            final spec = new ModuleSpec()
            spec.name = data.name as String
            spec.version = data.version as String
            spec.description = data.description as String
            spec.keywords = data.keywords as List<String> ?: []
            spec.license = data.license as String
            spec.authors = data.authors as List<String> ?: []
            spec.maintainers = data.maintainers as List<String> ?: []
            spec.requires = data.requires as Map<String, String> ?: [:]
            spec.tools = data.tools as List<Map> ?: []

            final inputs = data.input
            if( inputs instanceof List )
                spec.inputs = inputs.collect { moduleParam(it) }
            else if( inputs != null )
                throw new Exception("invalid spec inputs")

            final outputs = data.output
            if( outputs instanceof List )
                spec.outputs = outputs.collect { moduleParam(it) }
            else if( outputs instanceof Map )
                spec.outputs = outputs.entrySet().collect { e -> moduleParam(e.value) }
            else if( outputs != null )
                throw new Exception("invalid spec outputs")

            final topics = data.topics
            if( topics instanceof List )
                spec.topics = topics.collect { moduleParam(it) }
            else if( topics instanceof Map )
                spec.topics = topics.entrySet().collect { e -> moduleParam(e.value) }
            else if( topics != null )
                throw new Exception("invalid spec topics")

            spec._passthrough = data.subMap(data.keySet() - MODULE_SPEC_FIELDS)

            return spec
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to parse module spec: ${path}", e)
        }
    }

    /**
     * Load a module input/output from the YAML representation.
     *
     * Supports multiple structures:
     * - Tuples/records: list of param specs
     * - New module spec format: {name: "...", type: "..."}
     * - Old nf-core format: {<name>: {type: "...", ...}}
     *
     * @param item
     */
    private static ModuleParam moduleParam(Object item) {
        if( item instanceof List ) {
            // tuple/record param
            return moduleParamList(item)
        }

        if( item instanceof Map ) {
            // paramSpec or nf-core param
            return moduleParamMap(item)
        }

        throw new Exception("Invalid module param spec")
    }

    private static ModuleParam moduleParamList(List item) {
        final components = item.collect { moduleParam(it) }
        if( components.size() == 1 )
            return components[0]
        return new ModuleParam(
            name: null,
            type: null,
            components: components
        )
    }

    private static ModuleParam moduleParamMap(Map item) {
        if( item.size() == 1 && item.values().first() instanceof Map ) {
            // Old nf-core format: {<name>: {type: "...", description: "..."}}
            final name = item.keySet().first()
            final value = item[name] as Map
            return new ModuleParam(
                name: name,
                type: value['type'] as String,
                description: value['description'] as String,
                _passthrough: value.subMap(value.keySet() - PARAM_SPEC_FIELDS)
            )
        }
        else {
            // New paramSpec format: {name: "...", type: "..."}
            return new ModuleParam(
                name: item['name'] as String,
                type: item['type'] as String,
                description: item['description'] as String,
                _passthrough: item.subMap(item.keySet() - PARAM_SPEC_FIELDS)
            )
        }
    }

}
