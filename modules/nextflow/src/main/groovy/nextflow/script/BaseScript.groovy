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

import java.lang.reflect.InvocationTargetException
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.script.dsl.ProcessDslV1
import nextflow.script.dsl.ProcessDslV2
import nextflow.secret.SecretsLoader

/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class BaseScript extends Script implements ExecutionContext {

    private Session session

    private ScriptMeta meta

    private boolean typingEnabled

    private ParamsDef paramsDef

    private WorkflowDef entryFlow

    private OutputDef outputDef

    BaseScript() {
        meta = ScriptMeta.register(this)
    }

    BaseScript(Binding binding) {
        super(binding)
        meta = ScriptMeta.register(this)
    }

    @Override
    ScriptBinding getBinding() {
        (ScriptBinding)super.getBinding()
    }

    Session getSession() {
        session
    }

    boolean isTypingEnabled() {
        return typingEnabled
    }

    /**
     * The entry workflow of this script, or null if it doesn't have one.
     */
    WorkflowDef getEntryFlow() {
        return entryFlow
    }

    /**
     * The declared params of this script, keyed by name.
     */
    Map<String,Param> getParamDeclarations() {
        return paramsDef
            ? paramsDef.getDeclarations()
            : Collections.<String,Param>emptyMap()
    }

    /**
     * Whether this script defines an output block.
     *
     * The declared outputs are not resolved, because the output directives
     * of a pipeline can refer to its params, which are not resolved until
     * the pipeline is executed.
     */
    boolean hasOutputs() {
        return outputDef != null
    }

    /**
     * The outputs declared in the output block of this script, keyed by
     * name, with their output directives resolved against the given params.
     *
     * @param params
     */
    Map<String,Map> getOutputDeclarations(ScriptBinding.ParamsMap params) {
        return outputDef
            ? outputDef.getDeclarations(params)
            : Collections.<String,Map>emptyMap()
    }

    /**
     * Holds the configuration object which will used to execution the user tasks
     */
    @Deprecated
    protected Map getConfig() {
        final msg = "The access of `config` object is deprecated"
        throw new DeprecationException(msg)
    }

    /**
     * Enable disable task 'echo' configuration property
     * @param value
     */
    @Deprecated
    protected void echo(boolean value = true) {
        final msg = "The use of `echo` method has been deprecated"
        throw new DeprecationException(msg)
    }

    private void setup() {
        binding.owner = this
        session = binding.getSession()

        binding.setVariable( 'baseDir', session.baseDir )
        binding.setVariable( 'projectDir', session.baseDir )
        binding.setVariable( 'workDir', session.workDir )
        binding.setVariable( 'workflow', session.workflowMetadata )
        binding.setVariable( 'nextflow', NextflowMeta.instance )
        binding.setVariable( 'launchDir', Paths.get('./').toRealPath() )
        binding.setVariable( 'moduleDir', meta.moduleDir )
        binding.setVariable( 'secrets', SecretsLoader.secretContext() )
    }

    /**
     * Enable static typing for the script.
     */
    protected void enableTyping() {
        log.warn1 "Static typing is a preview feature -- syntax and behavior may change in future releases"
        this.typingEnabled = true
    }

    /**
     * Define a params block.
     *
     * @param clazz
     * @param body
     */
    protected void params(Class clazz, Closure body) {
        if( entryFlow )
            throw new IllegalStateException("Workflow params definition must be defined before the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow params definition is not allowed within a workflow")

        this.paramsDef = new ParamsDef(this, clazz, body)
    }

    /**
     * Define an agent.
     *
     * Mirrors {@link #processV2(String, Closure)} — the lowered agent closure runs
     * against an {@link AgentBuilder} delegate that captures directives, inputs,
     * outputs and the prompt, then builds the populated {@link AgentDef}.
     * The agent executes via {@link AgentDef#run} (see the nf-agent plugin runner).
     *
     * @param name
     * @param body
     */
    protected void agent(String name, Closure<PromptDef> body) {
        log.warn1 "Agents are a preview feature -- syntax and behavior may change in future releases"
        final builder = new AgentBuilder(this, name)
        final cl = (Closure<PromptDef>) body.clone()
        cl.setDelegate(builder)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        final prompt = cl.call()
        final agent = builder.withPrompt(prompt).build()
        meta.addDefinition(agent)
    }

    /**
     * Define a legacy process.
     *
     * @param name
     * @param body
     */
    protected void process(String name, Closure<BodyDef> body) {
        final dsl = new ProcessDslV1(this, name)
        final cl = (Closure<BodyDef>)body.clone()
        cl.setDelegate(dsl)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        final taskBody = cl.call()
        final process = dsl.withBody(taskBody).build()
        meta.addDefinition(process)
    }

    /**
     * Define a typed process.
     *
     * @param name
     * @param body
     */
    protected void processV2(String name, Closure<BodyDef> body) {
        final dsl = new ProcessDslV2(this, name)
        final cl = (Closure<BodyDef>)body.clone()
        cl.setDelegate(dsl)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        final taskBody = cl.call()
        final process = dsl.withBody(taskBody).build()
        meta.addDefinition(process)
    }

    /**
     * Define an entry workflow.
     *
     * @param workflowBody
     */
    protected void workflow(Closure<BodyDef> workflowBody) {
        // launch the execution
        final workflow = new WorkflowDef(this, workflowBody)
        // capture the main (unnamed) workflow definition
        this.entryFlow = workflow
        // add it to the list of workflow definitions
        meta.addDefinition(workflow)
    }

    /**
     * Define a named workflow.
     *
     * @param name
     * @param workflowBody
     */
    protected void workflow(String name, Closure<BodyDef> workflowBody) {
        final workflow = new WorkflowDef(this,workflowBody,name)
        meta.addDefinition(workflow)
    }

    /**
     * Define an output block.
     *
     * @param body
     */
    protected void output(Closure body) {
        if( !entryFlow )
            throw new IllegalStateException("Workflow output definition must be defined after the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow output definition is not allowed within a workflow")

        this.outputDef = new OutputDef(this, body)
    }

    /**
     * Include definitions from another script.
     *
     * @param include
     */
    protected IncludeDef include( IncludeDef include ) {
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Include statement is not allowed within a workflow definition")
        return include.setSession(session)
    }

    /**
     * Define a custom type.
     *
     * @param type
     */
    protected void declareType(Class type) {
        meta.addDefinition(new TypeDef(type))
    }

    /**
     * Invokes custom methods in the task execution context
     *
     * @see nextflow.processor.TaskContext#invokeMethod(java.lang.String, java.lang.Object)
     * @see WorkflowBinding#invokeMethod(java.lang.String, java.lang.Object)
     *
     * @param name the name of the method to call
     * @param args the arguments to use for the method call
     * @return The result of the custom method execution
     */
    @Override
    Object invokeMethod(String name, Object args) {
        binding.invokeMethod(name, args)
    }

    private Object run0() {
        final result = runScript()
        if( meta.isModule() ) {
            return result
        }

        // a module defines a single process or named workflow, which is
        // executed directly -- there is no entry workflow to select
        if( binding.entryName && session.isModuleRun() )
            throw new AbortOperationException("Option `-entry` is not supported when a script is executed as a module")

        // if an `entryName` was specified via the command line, override the `entryFlow` to be executed
        if( binding.entryName && !(entryFlow=meta.getWorkflow(binding.entryName) ) ) {
            def msg = "Unknown workflow entry name: ${binding.entryName}"
            final allNames = meta.getWorkflowNames()
            final guess = allNames.closest(binding.entryName)
            if( guess )
                msg += " -- Did you mean?\n" + guess.collect { "  $it"}.join('\n')
            throw new IllegalArgumentException(msg)
        }

        if( !entryFlow ) {
            // a process or named workflow can be executed directly only
            // when the script was launched via `nextflow module run`
            final moduleRun = session.isModuleRun()
            if( moduleRun && meta.hasExecutableWorkflows() ) {
                // Execute a single named workflow directly
                final handler = new WorkflowEntryHandler(this, session, meta)
                this.entryFlow = handler.createEntryWorkflow()
            }
            else if( moduleRun && meta.hasExecutableProcesses() ) {
                // Execute a single process directly
                final handler = new ProcessEntryHandler(this, session, meta)
                this.entryFlow = handler.createEntryWorkflow()
            }
            else if( meta.getLocalProcessNames() || meta.getLocalWorkflowNames() ) {
                throw new AbortOperationException("No entry workflow specified -- script must define an entry workflow, a single process or named workflow, or be a code snippet")
            }
            else {
                // NOTE: remove after v1 parser is removed
                return result
            }
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        // the entry workflow receives the resolved params as an input, so that
        // `params` refers to a single execution of the pipeline
        final params = paramsDef?.apply(session)
        if( params != null )
            entryFlow.withParams(params)
        final ret = entryFlow.invoke_a(BaseScriptConsts.EMPTY_ARGS)
        if( outputDef )
            outputDef.apply(session, params)
        session.notifyAfterWorkflowExecution()
        return ret
    }

    Object run() {
        setup()
        ExecutionStack.push(this)
        try {
            run0()
        }
        catch( InvocationTargetException e ) {
            // provide the exception cause which is more informative than InvocationTargetException
            Throwable target = e
            do target = target.cause
            while ( target instanceof InvocationTargetException )
            throw target
        }
        finally {
            ExecutionStack.pop()
        }
    }

    protected abstract Object runScript()

    @Override
    void print(Object object) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.print(object)
    }

    @Override
    void println() {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info("")
        else
            super.println()
    }

    @Override
    void println(Object object) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.println(object)
    }

    @Override
    void printf(String msg, Object arg) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(String.format(msg, arg))
        else
            super.printf(msg, arg)
    }

    @Override
    void printf(String msg, Object[] args) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(String.format(msg, args))
        else
            super.printf(msg, args)
    }

}
