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

    private ProcessFactory processFactory

    private ScriptMeta meta

    private WorkflowDef entryFlow

    private OutputDef publisher

    @Lazy InputStream stdin = { System.in }()

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
        processFactory = session.newProcessFactory(this)

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
     * Define a params block.
     *
     * @param body
     */
    protected void params(Closure body) {
        if( entryFlow )
            throw new IllegalStateException("Workflow params definition must be defined before the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow params definition is not allowed within a workflow")

        final dsl = new ParamsDsl()
        final cl = (Closure)body.clone()
        cl.setDelegate(dsl)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.call()

        dsl.apply(session)
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
     * @param closure
     */
    protected void output(Closure closure) {
        if( !NF.outputDefinitionEnabled )
            throw new IllegalStateException("Workflow output definition requires the `nextflow.preview.output` feature flag")
        if( !entryFlow )
            throw new IllegalStateException("Workflow output definition must be defined after the entry workflow")
        if( ExecutionStack.withinWorkflow() )
            throw new IllegalStateException("Workflow output definition is not allowed within a workflow")

        publisher = new OutputDef(closure)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(ExecutionStack.withinWorkflow())
            throw new IllegalStateException("Include statement is not allowed within a workflow definition")
        include .setSession(session)
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
            if( meta.getLocalWorkflowNames() )
                throw new AbortOperationException("No entry workflow specified")
            // Check if we have standalone processes that can be executed automatically
            if( meta.hasExecutableProcesses() ) {
                // Create a workflow to execute the process (single process or first of multiple)
                final handler = new ProcessEntryHandler(this, session, meta)
                entryFlow = handler.createAutoProcessEntry()
            }
            else {
                return result
            }
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        final ret = entryFlow.invoke_a(BaseScriptConsts.EMPTY_ARGS)
        if( publisher )
            publisher.apply(session)
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
