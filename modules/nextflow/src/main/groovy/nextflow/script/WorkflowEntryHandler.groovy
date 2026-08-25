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

package nextflow.script

import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type
import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovy.yaml.YamlSlurper
import nextflow.Session
import nextflow.dataflow.ChannelNamespace
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHelper
import nextflow.script.dsl.Types
import nextflow.script.types.Channel
import nextflow.script.types.Value
import nextflow.splitter.CsvSplitter
import nextflow.util.TypeHelper

/**
 * Helper class for named workflow execution.
 *
 * This feature enables direct execution of a named workflow without
 * an explicit entry workflow:
 * - Scripts with a single named workflow run it automatically:
 *   {@code nextflow run script.nf --param value}
 * - Command-line parameters are mapped directly to workflow inputs ({@code take:})
 * - Inputs of collection type are loaded from samplesheet files (CSV, JSON, YAML)
 * - Non-collection inputs are passed through as values
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class WorkflowEntryHandler {

    private final BaseScript script
    private final Session session
    private final WorkflowDef workflowDef

    WorkflowEntryHandler(BaseScript script, Session session, ScriptMeta meta) {
        this.script = script
        this.session = session

        final workflowNames = meta.getLocalWorkflowNames()
        if( workflowNames.size() != 1 )
            throw new IllegalStateException("Direct execution of named workflows is only supported for scripts with exactly one named workflow")

        if( !script.isTypingEnabled() )
            throw new IllegalStateException("Direct execution of named workflows is only supported when static typing is enabled")

        final workflowName = workflowNames.first()
        this.workflowDef = meta.getWorkflow(workflowName)
    }

    /**
     * Creates an entry workflow that calls the selected named workflow.
     *
     * Parameters are automatically mapped to workflow inputs, with
     * collection-typed inputs loaded from samplesheet files.
     *
     * Workflow emits are published as pipeline outputs, without creating
     * an output directory.
     */
    WorkflowDef createEntryWorkflow() {
        final workflowName = workflowDef.name
        final entryBody = { ->
            final entryExecutionClosure = { ->
                // Map parameters to workflow inputs
                final inputs = getWorkflowArguments(workflowDef)
                // Execute the named workflow
                final output = workflowDef.run(inputs as Object[]) as ChannelOut
                // Publish workflow emits as pipeline outputs
                assignOutputs(output)
                publishOutputs()
                return output
            }
            final sourceCode = "    // Auto-generated workflow entry\n    ${workflowName}(...)"
            return new BodyDef(entryExecutionClosure, sourceCode, 'workflow')
        }
        return new WorkflowDef(script, entryBody)
    }

    private void assignOutputs(ChannelOut output) {
        final outputNames = workflowDef.getDeclaredOutputs()
        final dsl = script.getBinding()
        if( output.size() == 1 && outputNames.size() == 1 ) {
            dsl._publish_(outputNames.first(), output[0])
        }
        else {
            for( final name : outputNames )
                dsl._publish_(name, output.getProperty(name))
        }
    }

    private void publishOutputs() {
        final outputNames = workflowDef.getDeclaredOutputs()
        final dsl = new OutputDsl()
        for( final name : outputNames )
            dsl.declare(name, { -> })
        // disable the output directory -- report output files by
        // their work directory path instead of publishing them
        session.outputDir = null
        dsl.apply(session)
    }

    /**
     * Resolves the workflow input arguments from the current session params.
     *
     * Each declared input ({@code take:} parameter) of the named workflow becomes
     * a pipeline parameter of the same name, and the declared type determines how
     * the parameter value is interpreted. This is the same mapping performed by
     * the {@code params} block for an entry workflow.
     *
     * @param workflowDef
     */
    protected List getWorkflowArguments(WorkflowDef workflowDef) {
        final inputs = workflowDef.getDeclaredInputs()
        final inputNames = inputs*.name
        final cliParams = session.cliParams ?: [:]
        final configParams = session.configParams ?: [:]

        for( final name : cliParams.keySet() ) {
            if( name !in inputNames && !configParams.containsKey(name) )
                throw new ScriptRuntimeException("Parameter `${name}` was specified on the command line but is not an input of workflow `${workflowDef.name}`")
        }

        final arguments = []
        for( final decl : inputs ) {
            final name = decl.name
            if( cliParams.containsKey(name) )
                arguments.add(resolveInput(decl, cliParams.get(name), true))
            else if( configParams.containsKey(name) )
                arguments.add(resolveInput(decl, configParams.get(name), false))
            else if( decl.optional )
                arguments.add(null)
            else
                throw new ScriptRuntimeException("Workflow `${workflowDef.name}` requires input `${name}` but no parameter `--${name}` was provided")
        }
        return arguments
    }

    /**
     * Resolves a single workflow input value against its declared type.
     *
     * A {@code Channel<E>} input is loaded from a samplesheet file, with each
     * record converted to the element type. A {@code Value<V>} input is
     * converted to {@code V} and wrapped in a value channel. Any other input is
     * converted directly to the declared type.
     *
     * @param decl    the declared input
     * @param value   the raw param value
     * @param fromCli whether the value came from the command line (and is
     *                therefore a string that may need to be parsed)
     */
    protected Object resolveInput(Param decl, Object value, boolean fromCli) {
        if( value == null )
            return null

        // an untyped input is passed through as-is
        if( decl.type == null )
            return value

        final rawType = TypeHelper.getRawType(decl.type)

        if( rawType == Channel )
            return ChannelNamespace.fromList(loadChannelInput(decl, value))

        if( rawType == Value )
            return ChannelNamespace.value(resolveScalarInput(elementDecl(decl), value, fromCli))

        return resolveScalarInput(decl, value, fromCli)
    }

    private Object resolveScalarInput(Param decl, Object value, boolean fromCli) {
        final result = fromCli
            ? ParamsHelper.resolveFromCli(decl, value)
            : ParamsHelper.resolveFromCode(decl, value)

        final expectedType = TypeHelper.getRawType(decl.type)
        final actualType = result?.getClass()
        if( actualType != null && !ParamsHelper.isAssignableFrom(expectedType, actualType) )
            throw new ScriptRuntimeException("Workflow input `${decl.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${result} [${Types.getName(actualType)}]")

        return result
    }

    /**
     * Loads a channel input from a samplesheet file, converting each record
     * to the declared element type.
     *
     * @param decl
     * @param value
     */
    private List loadChannelInput(Param decl, Object value) {
        if( value !instanceof CharSequence && value !instanceof Path )
            throw new ScriptRuntimeException("Workflow input `${decl.name}` with type ${Types.getName(decl.type)} should be a samplesheet file, but received: ${value} [${Types.getName(value.getClass())}]")

        final path = value instanceof Path
            ? (Path)value
            : FileHelper.asPath(value.toString())
        final elementType = elementDecl(decl).type
        return loadFromFile(decl.name, path).collect { el -> TypeHelper.asType(el, elementType) }
    }

    /**
     * Returns the declared input for the element type of a parameterized
     * input type, e.g. {@code Sample} for {@code Channel<Sample>}.
     *
     * @param decl
     */
    private static Param elementDecl(Param decl) {
        final elementType = decl.type instanceof ParameterizedType
            ? ((ParameterizedType)decl.type).getActualTypeArguments()[0]
            : (Type)Object
        return new Param(decl.name, elementType, decl.optional, null)
    }

    /**
     * Loads the contents of a samplesheet file as a list of records.
     *
     * Supported formats:
     * - CSV: header row required, comma-separated
     * - JSON: must be a top-level array
     * - YAML / YML: must be a top-level sequence
     *
     * @param name the input name (for error messages)
     * @param file the samplesheet file to load
     * @return a list of raw records (maps)
     */
    protected List loadFromFile(String name, Path file) {
        final ext = file.getExtension()
        final value = switch( ext ) {
            case 'csv'         -> new CsvSplitter().options(header: true, sep: ',').target(file).list()
            case 'json'        -> new JsonSlurper().parse(file)
            case 'yaml', 'yml' -> new YamlSlurper().parse(file)
            default -> throw new ScriptRuntimeException("Unrecognized file format '${ext}' for input file '${file}' for workflow input `${name}` -- should be CSV, JSON, or YAML")
        }
        if( value !instanceof List )
            throw new ScriptRuntimeException("Input file '${file}' for workflow input `${name}` must contain a list of records, but got: ${value.class.simpleName}")
        return (List)value
    }

}
