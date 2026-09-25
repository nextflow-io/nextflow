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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ScriptRuntimeException

/**
 * Helper class for named workflow execution.
 *
 * A script that defines a single named workflow and no processes can be
 * executed directly, without an explicit entry workflow:
 * {@code nextflow module run script.nf --param value}
 *
 * Each input ({@code take:}) becomes a pipeline parameter of the same name,
 * and its declared type determines how the param value is interpreted -- a
 * {@code Channel<E>} input is loaded from a samplesheet file (CSV, JSON,
 * YAML), while any other input is converted to the declared type. The
 * workflow emits are printed to standard output, without publishing them
 * to an output directory.
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

        final workflowName = workflowNames.first()
        if( !script.isTypingEnabled() )
            throw new ScriptRuntimeException("Workflow `${workflowName}` cannot be executed directly because it is not typed -- static typing is required to map pipeline parameters to workflow inputs")

        this.workflowDef = meta.getWorkflow(workflowName)

        // every input must be typed, because the declared type determines
        // how the corresponding param value is interpreted
        final untyped = workflowDef.getDeclaredInputs().findAll { decl -> decl.type == null }*.name
        if( untyped )
            throw new ScriptRuntimeException("Workflow `${workflowName}` cannot be executed directly because the following inputs are not typed: ${untyped.join(', ')}")
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
                assignOutputs((WorkflowBinding)(Object)getDelegate(), output)
                return output
            }
            final sourceCode = "    // Auto-generated workflow entry\n    ${workflowName}(...)"
            return new BodyDef(entryExecutionClosure, sourceCode, 'workflow')
        }
        return new WorkflowDef(script, entryBody)
    }

    private void assignOutputs(WorkflowBinding dsl, ChannelOut output) {
        final outputNames = workflowDef.getDeclaredOutputs()
        if( output.size() == 1 && outputNames.size() == 1 ) {
            dsl._publish_(outputNames.first(), output[0])
        }
        else {
            for( final name : outputNames )
                dsl._publish_(name, output.getProperty(name))
        }
    }

    /**
     * Creates an output definition that declares each output of the
     * entry workflow, without creating an output directory.
     */
    OutputDef createOutputDef() {
        final outputNames = workflowDef.getDeclaredOutputs()
        // disable the output directory -- report output files by
        // their work directory path instead of publishing them
        session.outputDir = null
        return new OutputDef({ ->
            final dsl = (OutputDsl)(Object)getDelegate()
            for( final name : outputNames )
                dsl.declare(name, { -> })
        })
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
            final value =
                cliParams.containsKey(name) ? ParamsHelper.resolveParam(decl, cliParams.get(name), true) :
                configParams.containsKey(name) ? ParamsHelper.resolveParam(decl, configParams.get(name), false) :
                ParamsHelper.resolveDefault(decl)

            if( value == null && !decl.optional ) {
                throw new ScriptRuntimeException("Parameter `--${name}` is required but no value was provided")
            }

            arguments.add(value)
        }
        return arguments
    }

}
