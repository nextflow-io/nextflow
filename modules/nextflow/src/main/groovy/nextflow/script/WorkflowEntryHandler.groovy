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

import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovy.yaml.YamlSlurper
import nextflow.Session
import nextflow.dataflow.ChannelNamespace
import nextflow.exception.ScriptRuntimeException
import nextflow.file.FileHelper
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
                final inputs = getWorkflowArguments(workflowDef, session.params)
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

    // TODO: disable output directory so that workflow output is produced without
    // actually copying files to output directory
    private void publishOutputs() {
        final outputNames = workflowDef.getDeclaredOutputs()
        final dsl = new OutputDsl()
        for( final name : outputNames )
            dsl.declare(name, { -> })
        dsl.apply(session)
    }

    /**
     * Resolves the workflow input arguments from the current session params.
     *
     * For each declared input ({@code take:} parameter) of the named workflow:
     * - If the input is typed as {@code Channel<T>}, the param value is resolved
     *   as a samplesheet path or collection and loaded into a channel
     * - If the param value is a collection, it is loaded into a channel via
     *   {@code channel.fromList()}
     * - If the param value is a string path to a samplesheet (CSV, JSON, YAML),
     *   the file is loaded and its contents are emitted as a channel
     * - Otherwise the value is passed directly as a workflow variable
     *
     * @param workflowDef
     * @param params
     */
    private List getWorkflowArguments(WorkflowDef workflowDef, Map params) {
        final inputs = workflowDef.getDeclaredInputs()
        final inputTypes = workflowDef.getDeclaredInputTypes()

        log.debug "Getting input arguments for workflow: ${workflowDef.name}"
        log.debug "Session params: ${params}"

        final arguments = []
        for( final name : inputs ) {
            final value = params.get(name)
            if( value == null && !params.containsKey(name) )
                throw new ScriptRuntimeException("Workflow `${workflowDef.name}` requires input `${name}` but no parameter `--${name}` was provided")

            final type = inputTypes.get(name)
            arguments.add(resolveInput(name, type, value))
        }

        log.debug "Final input arguments: ${arguments}"
        return arguments
    }

    /**
     * Resolves a single workflow input value.
     *
     * When the declared type is {@code Channel} (or a subtype), the value is
     * always resolved to a channel — either by loading a samplesheet file or by
     * wrapping an existing collection with {@code channel.fromList()}.
     *
     * For other (or unknown) types the value is passed through as-is, unless it
     * happens to be a collection or a samplesheet path, in which case a channel
     * is also created (heuristic fallback for untyped workflows).
     *
     * @param name  the input name (for error messages)
     * @param type  the declared type of the input, or {@code null} if untyped
     * @param value the raw param value
     * @return a {@code ChannelImpl} if the value resolves to a channel, or the
     *         raw value for scalar inputs
     */
    protected Object resolveInput(String name, Class type, Object value) {
        if( value == null )
            return value

        // TODO: need to generate __Params class to preserve input types
        if( type == Channel ) {
            if( value instanceof Collection ) {
                return ChannelNamespace.fromList((Collection)value)
            }
            if( value instanceof String ) {
                final path = FileHelper.asPath(value)
                return ChannelNamespace.fromList(loadFromFile(name, path))
            }
            throw new ScriptRuntimeException("Workflow input `${name}` expects a Channel but received: ${value} [${value.class.simpleName}]")
        }

        if( type == Value ) {
            return ChannelNamespace.value(value)
        }

        if( value !instanceof String )
            return TypeHelper.asType(value, type)

        final str = (String) value

        if( type == Boolean ) {
            if( str.toLowerCase() == 'true' ) return Boolean.TRUE
            if( str.toLowerCase() == 'false' ) return Boolean.FALSE
        }

        if( type == Integer || type == Float ) {
            if( str.isInteger() ) return str.toInteger()
            if( str.isLong() ) return str.toLong()
        }

        if( type == Float ) {
            if( str.isFloat() ) return str.toFloat()
            if( str.isDouble() ) return str.toDouble()
        }

        if( type == Path ) {
            return TypeHelper.asPathType(str)
        }

        return value
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
