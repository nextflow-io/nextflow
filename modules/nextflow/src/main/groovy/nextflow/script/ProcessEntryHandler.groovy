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
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.Nextflow
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpecFactory
import nextflow.module.ModuleStorage
import nextflow.script.dsl.Types
import nextflow.script.params.DefaultInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.InParam
import nextflow.script.params.StdInParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessTupleInput
import nextflow.script.types.Record
import nextflow.util.RecordMap

/**
 * Helper class for process entry execution feature.
 *
 * A script that defines a single process and no workflows can be executed
 * directly, without an explicit entry workflow:
 * {@code nextflow module run script.nf --param value}
 *
 * Command-line parameters are mapped to the process inputs, and the process
 * outputs become the workflow outputs. Supports the following process input
 * qualifiers: val, path, tuple, each
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessEntryHandler {

    private final BaseScript script
    private final Session session
    private final ProcessDef processDef

    ProcessEntryHandler(BaseScript script, Session session, ScriptMeta meta) {
        this.script = script
        this.session = session

        final processNames = meta.getLocalProcessNames()
        if( processNames.size() != 1 )
            throw new IllegalStateException("Direct execution of processes is only supported for scripts with exactly one process")

        final processName = processNames.first()
        this.processDef = meta.getProcess(processName)
    }

    /**
     * Creates a workflow to execute the specified process.
     *
     * Parameters are automatically mapped to process inputs, and
     * process outputs are mapped to workflow outputs.
     */
    WorkflowDef createEntryWorkflow() {
        final processName = processDef.name
        final workflowBody = { ->
            // Create the workflow execution logic
            final workflowExecutionClosure = { ->
                // Map parameters to process inputs
                final inputArgs = getProcessArguments(processDef)
                // Execute the process
                final output = processDef.run(inputArgs as Object[]) as ChannelOut
                // Publish process outputs as workflow outputs
                assignOutputs(output)
                publishOutputs()
                return output
            }

            // Create the body definition with execution logic
            final sourceCode = "    // Auto-generated process entry\n    ${processName}(...)"
            return new BodyDef(workflowExecutionClosure, sourceCode, 'workflow')
        }

        return new WorkflowDef(script, workflowBody)
    }

    private void assignOutputs(ChannelOut output) {
        final config = processDef.getProcessConfig()
        final outputNames = getProcessOutputs(config)
        final dsl = script.getBinding()
        if( output.size() != outputNames.size() )
            log.warn("Process ${processDef.name} is missing emit names for one or more outputs -- unnamed outputs will be omitted")
        if( output.size() == 1 && outputNames.size() == 1 ) {
            dsl._publish_(outputNames.first(), output[0])
        }
        else {
            for( final name : outputNames )
                dsl._publish_(name, output.getProperty(name))
        }
    }

    private void publishOutputs() {
        final config = processDef.getProcessConfig()
        final outputNames = getProcessOutputs(config)
        final dsl = new OutputDsl()
        for( final name : outputNames )
            dsl.declare(name, { -> })
        // disable the output directory -- report output files by
        // their work directory path instead of publishing them
        session.outputDir = null
        dsl.apply(session)
    }

    private List<String> getProcessOutputs(ProcessConfig config) {
        return config instanceof ProcessConfigV1
            ? getProcessOutputsV1(config)
            : getProcessOutputsV2((ProcessConfigV2) config)
    }

    private List<String> getProcessOutputsV1(ProcessConfigV1 config) {
        final output = config.getOutputs()
        final names = output*.getChannelEmitName().findAll { it != null }
        if( output.isEmpty() )
            return Collections.emptyList()
        if( output.size() == 1 && names.isEmpty() )
            return List.of('$out')
        return names
    }

    private List<String> getProcessOutputsV2(ProcessConfigV2 config) {
        final output = config.getOutputs().getParams()
        final names = output*.getName()
        if( output.isEmpty() )
            return Collections.emptyList()
        if( output.size() == 1 && (names.isEmpty() || names.first() == '$out') )
            return List.of('$out')
        return names
    }

    /**
     * Gets the input arguments for a process by mapping the params given on
     * the command line and in the config to the declared process inputs.
     *
     * <p>Uses the same binding as the reusable static
     * {@link #getProcessArguments(ProcessDef, Map, ModuleSpec)}, passing the sibling
     * {@code meta.yml} resolved from the running script path so the existing
     * {@code module run} behavior (dot-params, type coercion from {@code meta.yml},
     * tuple assembly) is preserved exactly.
     *
     * <p>A typed (V2) input is resolved from the command-line params, then the
     * remaining params, in the same way as a pipeline param or a workflow
     * input, since a command-line value is a string that must be parsed
     * whereas a config or script value is already structured.
     *
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    protected List getProcessArguments(ProcessDef processDef) {
        final scriptPath = script?.getBinding()?.getScriptPath()
        return bindProcessArguments(processDef, getModuleSpecInputTypes(scriptPath), session.params ?: [:], session.cliParams ?: [:])
    }

    /**
     * Maps {@code params} onto the declared inputs of {@code processDef}, returning one element
     * per input channel (a tuple input becomes a {@code List} of its component values, e.g.
     * {@code [[id:'s1'], file(reads)]}). This is the reusable, instance-free form of the
     * {@code module run} param→channel binding: dot-notation params are folded into nested maps,
     * legacy (V1) inputs are coerced using the input TYPES declared in the given module spec,
     * and typed (V2) inputs are coerced from their declared
     * {@link nextflow.script.params.v2.ProcessInput} type. Every value is treated as
     * command-line text.
     *
     * @param processDef the process whose inputs are bound
     * @param params     the (possibly dotted) param map to bind by input name
     * @param spec       the module spec providing input types for the legacy (V1) path; may be
     *                   {@code null}, in which case an empty type map is used
     * @return list of values to pass to {@link ProcessDef#run}, one per input channel
     */
    static List getProcessArguments(ProcessDef processDef, Map params, ModuleSpec spec) {
        return bindProcessArguments(processDef, spec != null ? moduleSpecInputTypes(spec) : Collections.<String,Class>emptyMap(), params, params)
    }

    /**
     * @param processDef the process whose inputs are bound
     * @param paramTypes the input types declared in the module spec, for the legacy (V1) path
     * @param params     all param values, as structured by the config and the script
     * @param cliParams  the param values given on the command line or in a params file,
     *                   as raw strings; the same map as {@code params} when every value
     *                   should be treated as such
     */
    private static List bindProcessArguments(ProcessDef processDef, Map<String,Class> paramTypes, Map params, Map cliParams) {
        try {
            log.debug "Getting input arguments for process: ${processDef.name}"
            log.debug "Session params: ${params}"

            final config = processDef.getProcessConfig()
            final inputArgs = config instanceof ProcessConfigV1
                ? getProcessArgumentsV1(config, params, paramTypes)
                : getProcessArgumentsV2((ProcessConfigV2) config, params, cliParams)

            log.debug "Final input arguments: ${inputArgs}"
            return inputArgs

        } catch (Exception e) {
            log.error "Failed to get input arguments for process ${processDef.name}: ${e.message}"
            throw e
        }
    }

    private static List getProcessArgumentsV1(ProcessConfigV1 config, Map params, Map<String,Class> paramTypes) {
        final declaredInputs = config.getInputs()

        if( declaredInputs.isEmpty() ) {
            return []
        }

        // Map declared inputs to command-line arguments
        List arguments = []
        for( final param : declaredInputs ) {
            // Skip the synthetic `$` control input that a process gains once it has been
            // `run()` (DefaultInParam): it is a termination-control channel, never a
            // user-supplied value. It is absent in the typical `module run` path (which binds
            // BEFORE run) and present when binding a process that was pre-wired/run earlier
            // (e.g. the agent tool bridge) - skipping it makes both paths produce the same args.
            if( param instanceof DefaultInParam )
                continue
            if( param instanceof TupleInParam ) {
                List tupleElements = []
                for( final innerParam : param.inner ) {
                    final value = getValueForInputV1(innerParam, params, paramTypes)
                    tupleElements.add(value)
                }
                arguments.add(tupleElements)
            }
            else {
                final value = getValueForInputV1(param, params, paramTypes)
                arguments.add(value)
            }
        }

        return arguments
    }

    /**
     * Load mapping of input types from the module spec if available. Returns empty map
     * if module spec is absent or unreadable.
     */
    private static Map<String, Class> getModuleSpecInputTypes(Path scriptPath) {
        try {
            final specPath = scriptPath?.resolveSibling(ModuleStorage.MODULE_MANIFEST_FILE)
            final spec = ModuleSpecFactory.fromYaml(specPath)
            return moduleSpecInputTypes(spec)
        }
        catch( Exception e ) {
            log.debug "Could not load input types from module spec: ${e.message}"
        }
        return Collections.emptyMap()
    }

    private static Map<String, Class> moduleSpecInputTypes(ModuleSpec spec) {
        final Map<String, Class> result = [:]
        for( final input : spec.inputs ) {
            if( input.isTuple() ) {
                for( final el : input.components )
                    result.put(el.name, inputType(el.type))
            }
            else {
                result.put(input.name, inputType(input.type))
            }
        }
        return result
    }

    private static Class inputType(String type) {
        return switch( type ) {
            case 'boolean' -> Boolean
            case 'file' -> Path
            case 'directory' -> Path
            case 'float' -> Number
            case 'integer' -> Integer
            case 'map' -> Map
            case 'string' -> String
            default -> null
        }
    }

    /**
     * Gets the appropriate value for a legacy process input.
     *
     * @param decl Input declaration
     * @param namedArgs Map of command-line arguments
     * @param paramTypes Map of input types from module spec
     * @return Properly typed value for the input
     */
    private static Object getValueForInputV1(InParam decl, Map namedArgs, Map<String,Class> paramTypes) {
        final name = decl.getName()
        final type = paramTypes.get(name)
        final value = namedArgs.get(name)

        // File/path inputs: an ABSENT value means "not provided". nf-core path inputs are
        // optional by convention and default to an empty list (the process script handles the
        // empty case). An empty value is NOT a stand-in for an absent one: as with most CLI
        // tools, an optional path input is skipped by supplying nothing at all, not by
        // supplying the option with an empty value. An empty value therefore falls through to
        // `file('')`, which fails loudly.
        if( decl instanceof FileInParam ) {
            if( value == null ) {
                log.warn "Path input '--${name}' not provided, defaulting to empty list"
                return []
            }
            return parseFileInput(value.toString())
        }

        // non-file inputs: a missing value is a hard error (required)
        if( value == null )
            throw new IllegalArgumentException("Parameter `--${name}` is required but no value was provided")

        // handle env, stdin inputs
        switch( decl ) {
            case EnvInParam:
                throw new IllegalArgumentException("Process `env` input qualifier is not supported by implicit process entry")

            case StdInParam:
                throw new IllegalArgumentException("Process `stdin` input qualifier is not supported by implicit process entry")
        }

        // handle val inputs
        if( value !instanceof String ) {
            if( type != null && !type.isAssignableFrom(value.class) )
                throw new IllegalArgumentException("Parameter '--${name}' expects a ${type.simpleName} but received: ${value} [${value.class.simpleName}]")
            return value
        }

        final str = (String) value

        if( type == Boolean ) {
            if( str.toLowerCase() == 'true' ) return Boolean.TRUE
            if( str.toLowerCase() == 'false' ) return Boolean.FALSE
            throw new IllegalArgumentException("Parameter '--${name}' expects a boolean (true/false) but received: '${str}'")
        }

        if( type == Integer ) {
            if( str.isInteger() ) return str.toInteger()
            if( str.isLong() ) return str.toLong()
            throw new IllegalArgumentException("Parameter '--${name}' expects an integer but received: '${str}'")
        }

        if( type == Number ) {
            if( str.isFloat() ) return str.toFloat()
            if( str.isDouble() ) return str.toDouble()
            throw new IllegalArgumentException("Parameter '--${name}' expects a floating-point number but received: '${str}'")
        }

        if( type == Map ) {
            throw new IllegalArgumentException(
                "Parameter '--${name}' expects a map but received: '${str}'. " +
                "Use dot notation (e.g. --${name}.key=value) to supply map properties.")
        }

        return str
    }

    private static List getProcessArgumentsV2(ProcessConfigV2 config, Map params, Map cliParams) {
        final declaredInputs = config.getInputs().getParams()

        if( declaredInputs.isEmpty() ) {
            return []
        }

        // Map declared inputs to command-line arguments
        List arguments = []
        for( final param : declaredInputs ) {
            if( param instanceof ProcessTupleInput && param.getType() == Record.class ) {
                final Map<String,Object> recordFields = [:]
                for( final innerParam : param.getComponents() ) {
                    final value = getValueForInputV2(innerParam, params, cliParams)
                    recordFields.put(innerParam.getName(), value)
                }
                arguments.add(new RecordMap(recordFields))
            }
            else if( param instanceof ProcessTupleInput ) {
                final List tupleElements = []
                for( final innerParam : param.getComponents() ) {
                    final value = getValueForInputV2(innerParam, params, cliParams)
                    tupleElements.add(value)
                }
                arguments.add(tupleElements)
            }
            else {
                final value = getValueForInputV2(param, params, cliParams)
                arguments.add(value)
            }
        }

        return arguments
    }

    /**
     * Gets the appropriate value for a typed process input.
     *
     * A typed input is resolved in the same way as a pipeline parameter,
     * so that a process and a workflow accept the same command line.
     *
     * @param input Input declaration
     * @param params All param values, as structured by the config and the script
     * @param cliParams The param values given on the command line or in a params file
     * @return Properly typed value for the input
     */
    private static Object getValueForInputV2(ProcessInput input, Map params, Map cliParams) {
        final name = input.getName()
        final type = input.getType()
        final decl = new Param(name, type, input.isOptional(), null)

        // a command-line value is a string that must be parsed, whereas a
        // config or script value is already structured
        final result =
            cliParams.containsKey(name) ? ParamsHelper.resolveFromCli(decl, cliParams.get(name)) :
            params.containsKey(name) ? ParamsHelper.resolveFromCode(decl, params.get(name)) :
            null

        if( result == null ) {
            if( decl.isOptional() )
                return null
            throw new IllegalArgumentException("Parameter `--${name}` is required but no value was provided")
        }

        // report a value that could not be converted
        final actualType = result?.getClass()
        if( type != null && actualType != null && !ParamsHelper.isAssignableFrom(type, actualType) )
            throw new IllegalArgumentException("Parameter `--${name}` with type ${Types.getName(type)} cannot be assigned to ${result} [${Types.getName(actualType)}]")

        return result
    }

    /**
     * Parses file input handling comma-separated values.
     * If the input contains commas, splits and returns a list of files.
     * Otherwise returns a single file object.
     *
     * @param fileInput String representation of file path(s)
     * @return Single file or list of files
     */
    protected static Object parseFileInput(String fileInput) {
        if( fileInput.contains(',') ) {
            // Split by comma, trim whitespace, and convert each to a file
            return fileInput.tokenize(',')
                .collect { it.trim() }
                .findAll { !it.isEmpty() }
                .collect { Nextflow.file(it) }
        } else {
            // Single file case - existing behavior
            return Nextflow.file(fileInput)
        }
    }
}
