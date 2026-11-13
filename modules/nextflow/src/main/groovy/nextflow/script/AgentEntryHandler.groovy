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
import nextflow.Session
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessTupleInput

/**
 * Helper class for agent entry execution feature.
 *
 * A script that defines a single agent and no workflows can be executed
 * directly, without an explicit entry workflow:
 * {@code nextflow module run script.nf --param value}
 *
 * Parameters are mapped to the agent inputs in the same way as a typed
 * process, and the agent outputs become the workflow outputs. Processes
 * may be defined alongside the agent, since they can be used as tools.
 */
@CompileStatic
class AgentEntryHandler {

    private final BaseScript script
    private final Session session
    private final AgentDef agentDef

    AgentEntryHandler(BaseScript script, Session session, ScriptMeta meta) {
        this.script = script
        this.session = session

        final agentNames = meta.getLocalAgentNames()
        if( agentNames.size() != 1 )
            throw new IllegalStateException("Direct execution of agents is only supported for scripts with exactly one agent")

        this.agentDef = meta.getAgent(agentNames.first())
    }

    /**
     * Creates a workflow to execute the agent.
     */
    WorkflowDef createEntryWorkflow() {
        final agentName = agentDef.name
        final workflowBody = { ->
            final workflowExecutionClosure = { ->
                final output = agentDef.run(getAgentArguments() as Object[]) as ChannelOut
                final dsl = (WorkflowBinding)(Object)getDelegate()
                dsl._publish_(agentDef.getOutput().name, output[0])
                return output
            }

            final sourceCode = "    // Auto-generated agent entry\n    ${agentName}(...)"
            return new BodyDef(workflowExecutionClosure, sourceCode, 'workflow')
        }

        return new WorkflowDef(script, workflowBody)
    }

    /**
     * Creates an output definition that declares the output of the
     * entry workflow, without creating an output directory.
     */
    OutputDef createOutputDef() {
        session.outputDir = null
        return new OutputDef({ ->
            final dsl = (OutputDsl)(Object)getDelegate()
            dsl.declare(agentDef.getOutput().name, { -> })
        })
    }

    /**
     * Gets the input arguments for the agent by mapping the params given
     * on the command line and in the config to the declared agent inputs.
     */
    protected List getAgentArguments() {
        final List<ProcessInput> inputs = agentDef.getInputs().collect { input ->
            input.components != null
                ? (ProcessInput) new ProcessTupleInput(input.components, input.type)
                : new ProcessInput(input.name, input.type, input.optional)
        }
        return ProcessEntryHandler.getProcessArgumentsV2(inputs, session.params ?: [:], session.cliParams ?: [:])
    }
}
