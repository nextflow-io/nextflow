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
import nextflow.module.ModuleStorage
import nextflow.script.params.EnvInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.InParam
import nextflow.script.params.StdInParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessTupleInput

/**
 * Helper class for process entry execution feature.
 *
 * This feature enables direct execution of Nextflow processes without explicit workflows:
 * - Single process scripts run automatically: `nextflow run script.nf --param value`
 * - Multi-process scripts run the first process automatically: `nextflow run script.nf --param value`
 * - Command-line parameters are mapped directly to process inputs
 * - Process outputs are mapped to workflow outputs
 * - Supports the following process input qualifiers: val, path, tuple, each
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
        if( processNames.isEmpty() )
            throw new IllegalStateException("No processes found for automatic execution")

        // Always pick the first process (whether single or multiple processes)
        final processName = processNames.first()
        this.processDef = meta.getProcess(processName)
    }

    /**
     * Creates a workflow to execute the specified process.
     *
     * Parameters are autoamtically mapped to process inputs, and
     * process outputs are mapped to workflow outputs.
     */
    WorkflowDef createEntryWorkflow() {
        final processName = processDef.name
        final workflowBody = { ->
            // Create the workflow execution logic
            final workflowExecutionClosure = { ->
                // Map parameters to process inputs
                final inputArgs = getProcessArguments(processDef, session.params)
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
            return List.of('result')
        return names
    }

    private List<String> getProcessOutputsV2(ProcessConfigV2 config) {
        final output = config.getOutputs().getParams()
        final names = output*.getName()
        if( output.isEmpty() )
            return Collections.emptyList()
        if( output.size() == 1 && (names.isEmpty() || names.first() == '$out') )
            return List.of('result')
        return names
    }

    /**
     * Gets the input arguments for a process by parsing input parameter structures
     * and mapping them from session.params, supporting dot notation for complex inputs.
     *
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    private List getProcessArguments(ProcessDef processDef, Map params) {
        try {
            log.debug "Getting input arguments for process: ${processDef.name}"
            log.debug "Session params: ${params}"

            final config = processDef.getProcessConfig()
            final inputArgs = config instanceof ProcessConfigV1
                ? getProcessArgumentsV1(config, params)
                : getProcessArgumentsV2((ProcessConfigV2) config, params)

            log.debug "Final input arguments: ${inputArgs}"
            return inputArgs

        } catch (Exception e) {
            log.error "Failed to get input arguments for process ${processDef.name}: ${e.message}"
            throw e
        }
    }

    /**
     * Parses complex parameters with dot notation support.
     * Converts flat parameters like --meta.id=1 --meta.name=test to nested maps.
     *
     * @param params Flat parameter map from session.params
     * @return Map with nested structures for complex parameters
     */
    private Map parseComplexParameters(Map params) {
        Map complexParams = [:]

        for( final entry : params.entrySet() ) {
            final parts = entry.key.toString().split('\\.')
            if( parts.length > 1 ) {
                Map current = complexParams
                for( int i = 0; i < parts.length - 1; i++ ) {
                    if( current[parts[i]] !instanceof Map ) {
                        current[parts[i]] = [:]
                    }
                    current = (Map) current[parts[i]]
                }
                current[parts[-1]] = entry.value
            } else {
                complexParams[entry.key] = entry.value
            }
        }

        return complexParams
    }

    private List getProcessArgumentsV1(ProcessConfigV1 config, Map params) {
        final declaredInputs = config.getInputs()

        if( declaredInputs.isEmpty() ) {
            return []
        }

        // Parse parameter values from session.params (handles dot notation)
        final paramValues = parseComplexParameters(params)

        // Load parameter types from module spec (if available)
        final paramTypes = getModuleSpecInputTypes()

        log.debug "Parameter values: ${paramValues}"

        // Map declared inputs to command-line arguments
        List arguments = []
        for( final param : declaredInputs ) {
            if( param instanceof TupleInParam ) {
                List tupleElements = []
                for( final innerParam : param.inner ) {
                    final value = getValueForInputV1(innerParam, paramValues, paramTypes)
                    tupleElements.add(value)
                }
                arguments.add(tupleElements)
            }
            else {
                final value = getValueForInputV1(param, paramValues, paramTypes)
                arguments.add(value)
            }
        }

        return arguments
    }

    /**
     * Load mapping of input types from the module spec if available. Returns null
     * if module spec is absent or unreadable.
     */
    private Map<String, Class> getModuleSpecInputTypes() {
        try {
            final scriptPath = script.getBinding().getScriptPath()
            final specPath = scriptPath?.resolveSibling(ModuleStorage.MODULE_MANIFEST_FILE)
            if( specPath )
                return ModuleSpec.loadInputTypes(specPath)
        }
        catch( Exception e ) {
            log.debug "Could not load input types from module spec: ${e.message}"
        }
        return Collections.emptyMap()
    }

    /**
     * Gets the appropriate value for a legacy process input.
     *
     * @param param Input declaration
     * @param namedArgs Map of command-line arguments
     * @param paramTypes Map of input types from module spec
     * @return Properly typed value for the input
     */
    private Object getValueForInputV1(InParam param, Map namedArgs, Map<String,Class> paramTypes) {
        final name = param.getName()
        final type = paramTypes.get(name)
        final value = namedArgs.get(name)

        if( value == null ) {
            throw new IllegalArgumentException("Missing required parameter: --${name}")
        }

        // handle file, path, env, stdin inputs
        switch( param ) {
            case FileInParam:
                return parseFileInput(value.toString())

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

    private List getProcessArgumentsV2(ProcessConfigV2 config, Map params) {
        final declaredInputs = config.getInputs().getParams()

        if( declaredInputs.isEmpty() ) {
            return []
        }

        // Parse complex parameters from session.params (handles dot notation)
        final paramValues = parseComplexParameters(params)
        log.debug "Parameter values: ${paramValues}"

        // Map declared inputs to command-line arguments
        List arguments = []
        for( final param : declaredInputs ) {
            if( param instanceof ProcessTupleInput ) {
                List tupleElements = []
                for( final innerParam : param.getComponents() ) {
                    final value = getValueForInputV2(innerParam, paramValues)
                    tupleElements.add(value)
                }
                arguments.add(tupleElements)
            }
            else {
                final value = getValueForInputV2(param, paramValues)
                arguments.add(value)
            }
        }

        return arguments
    }

    /**
     * Gets the appropriate value for a typed process input.
     *
     * @param param Input declaration
     * @param namedArgs Map of command-line arguments
     * @return Properly typed value for the input
     */
    private Object getValueForInputV2(ProcessInput param, Map namedArgs) {
        final name = param.getName()
        final type = param.getType()
        final value = namedArgs.get(name)

        if( value == null ) {
            throw new IllegalArgumentException("Missing required parameter: --${name}")
        }

        if( value !instanceof String )
            return value

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
            return Nextflow.file(str)
        }

        return value
    }

    /**
     * Parses file input handling comma-separated values.
     * If the input contains commas, splits and returns a list of files.
     * Otherwise returns a single file object.
     *
     * @param fileInput String representation of file path(s)
     * @return Single file object or List of file objects
     */
    protected Object parseFileInput(String fileInput) {
        if (fileInput.contains(',')) {
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
