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

package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.Nextflow
import nextflow.script.params.EnvInParam
import nextflow.script.params.FileInParam
import nextflow.script.params.InParam
import nextflow.script.params.StdInParam
import nextflow.script.params.TupleInParam

/**
 * Helper class for process entry execution feature.
 *
 * This feature enables direct execution of Nextflow processes without explicit workflows:
 * - Single process scripts run automatically: `nextflow run script.nf --param value`
 * - Multi-process scripts run the first process automatically: `nextflow run script.nf --param value`
 * - Command-line parameters are mapped directly to process inputs
 * - Supports the following process input qualifiers: val, path, tuple, each
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessEntryHandler {

    private final BaseScript script
    private final Session session
    private final ScriptMeta meta

    ProcessEntryHandler(BaseScript script, Session session, ScriptMeta meta) {
        this.script = script
        this.session = session
        this.meta = meta
    }

    /**
     * Creates a workflow to execute a standalone process automatically.
     * For single process scripts, executes that process.
     * For multi-process scripts, executes the first process.
     *
     * @return WorkflowDef that executes the process with parameter mapping
     */
    WorkflowDef createAutoProcessEntry() {
        final processNames = meta.getLocalProcessNames()
        if( processNames.isEmpty() ) {
            throw new IllegalStateException("No processes found for auto-execution")
        }

        // Always pick the first process (whether single or multiple processes)
        final processName = processNames.first()
        final processDef = meta.getProcess(processName)

        return createProcessEntry(processDef)
    }

    /**
     * Creates a workflow to execute the specified process with automatic parameter mapping.
     */
    private WorkflowDef createProcessEntry(ProcessDef processDef) {
        final processName = processDef.name

        final workflowBody = { ->
            // Create the workflow execution logic
            final workflowExecutionClosure = { ->
                // Get input parameter values and execute the process
                final inputArgs = getProcessInputArguments(processDef)
                final processResult = script.invokeMethod(processName, inputArgs as Object[])

                return processResult
            }

            // Create the body definition with execution logic
            final sourceCode = "    // Auto-generated process entry\n    ${processName}(...)"
            return new BodyDef(workflowExecutionClosure, sourceCode, 'workflow')
        }

        return new WorkflowDef(script, workflowBody)
    }

    /**
     * Gets the input arguments for a process by parsing input parameter structures
     * and mapping them from session.params, supporting dot notation for complex inputs.
     *
     * @param processDef The ProcessDef object containing the process definition
     * @return List of parameter values to pass to the process
     */
    private List getProcessInputArguments(ProcessDef processDef) {
        try {
            log.debug "Getting input arguments for process: ${processDef.name}"
            log.debug "Session params: ${session.params}"

            final declaredInputs = processDef.getProcessConfig().getInputs()

            if( declaredInputs.isEmpty() ) {
                return []
            }

            // Parse complex parameters from session.params (handles dot notation)
            final complexParams = parseComplexParameters(session.params)
            log.debug "Complex parameters: ${complexParams}"

            // Map declared inputs to command-line arguments
            List inputArgs = []
            for( final param : declaredInputs ) {
                if( param instanceof TupleInParam ) {
                    List tupleElements = []
                    for( final innerParam : param.inner ) {
                        final value = getValueForInput(innerParam, complexParams)
                        tupleElements.add(value)
                    }
                    inputArgs.add(tupleElements)
                }
                else {
                    final value = getValueForInput(param, complexParams)
                    inputArgs.add(value)
                }
            }

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

        params.each { key, value ->
            def parts = key.toString().split('\\.')
            if( parts.length > 1 ) {
                // Handle dot notation - build complex map
                def current = complexParams
                for( int i = 0; i < parts.length - 1; i++ ) {
                    if( !current.containsKey(parts[i]) ) {
                        current[parts[i]] = [:]
                    }
                    current = current[parts[i]]
                }
                current[parts[-1]] = value
            } else {
                // Simple parameter
                complexParams[key] = value
            }
        }

        return complexParams
    }

    /**
     * Gets the appropriate value for an input definition, handling type conversion.
     *
     * @param param Input declaration
     * @param namedArgs Map of command-line arguments
     * @return Properly typed value for the input
     */
    private Object getValueForInput(InParam param, Map namedArgs) {
        final paramName = param.getName()
        final paramValue = namedArgs.get(paramName)

        if( paramValue == null ) {
            throw new IllegalArgumentException("Missing required parameter: --${paramName}")
        }

        switch( param ) {
            case FileInParam:
                if( paramValue instanceof String ) {
                    return Nextflow.file(paramValue)
                }
                return paramValue

            case EnvInParam:
                throw new IllegalArgumentException("Process `env` input qualifier is not supported by implicit process entry")

            case StdInParam:
                throw new IllegalArgumentException("Process `stdin` input qualifier is not supported by implicit process entry")

            default:
                return paramValue
        }
    }
}
