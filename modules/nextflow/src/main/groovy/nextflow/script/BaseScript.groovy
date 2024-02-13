/*
 * Copyright 2013-2023, Seqera Labs
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

import groovy.transform.PackageScope
import nextflow.processor.TaskProcessor

import java.lang.reflect.InvocationTargetException
import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.exception.AbortOperationException
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script implements ExecutionContext {

    private TaskProcessor taskProcessor

    private Session session

    private ProcessFactory processFactory

    private ScriptMeta meta

    private WorkflowDef entryFlow

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

    /**
     * Access to the last *process* object -- only for testing purpose
     */
    @PackageScope
    TaskProcessor getTaskProcessor() { taskProcessor }

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
        binding.setVariable('launchDir', Paths.get('./').toRealPath())
        binding.setVariable('moduleDir', meta.moduleDir )
    }

    protected process( String name, Closure<BodyDef> body ) {
        if( NF.isDsl2() ) {
            def process = new ProcessDef(this,body,name)
            meta.addDefinition(process)
        }
        else {
            // legacy process definition an execution
            taskProcessor = processFactory.createProcessor(name, body)
            taskProcessor.run()
        }
    }

    /**
     * Workflow main entry point
     *
     * @param body The implementation body of the workflow
     * @return The result of workflow execution
     */
    protected workflow(Closure<BodyDef> workflowBody) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.enable.dsl=2` to allow the definition of workflow components")
        // launch the execution
        final workflow = new WorkflowDef(this, workflowBody)
        // capture the main (unnamed) workflow definition
        this.entryFlow = workflow
        // add it to the list of workflow definitions
        meta.addDefinition(workflow)
    }

    protected workflow(String name, Closure<BodyDef> workflowDef) {
        final workflow = new WorkflowDef(this,workflowDef,name)
        meta.addDefinition(workflow)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.enable.dsl=2` to import module files")
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
        if(NF.isDsl2())
            binding.invokeMethod(name, args)
        else
            super.invokeMethod(name, args)
    }

    private runDsl1() {
        session.notifyBeforeWorkflowExecution()
        final ret = runScript()
        session.notifyAfterWorkflowExecution()
        return ret
    }

    private runDsl2() {
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
                log.warn "No entry workflow specified"
            if( meta.getLocalProcessNames() ) {
                final msg = """\
                        =============================================================================
                        =                                WARNING                                    =
                        = You are running this script using DSL2 syntax, however it does not        = 
                        = contain any 'workflow' definition so there's nothing for Nextflow to run. =
                        =                                                                           =
                        = If this script was written using Nextflow DSL1 syntax, please add the     = 
                        = setting 'nextflow.enable.dsl=1' to the nextflow.config file or use the    =
                        = command-line option '-dsl1' when running the pipeline.                    =
                        =                                                                           =
                        = More details at this link: https://www.nextflow.io/docs/latest/dsl2.html  =
                        =============================================================================
                        """.stripIndent(true)
                throw new AbortOperationException(msg)
            }
            return result
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        final ret = entryFlow.invoke_a(BaseScriptConsts.EMPTY_ARGS)
        session.notifyAfterWorkflowExecution()
        return ret
    }

    Object run() {
        setup()
        ExecutionStack.push(this)
        try {
//            run0()
            NF.dsl2 ? runDsl2() : runDsl1()
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
            log.info(String.printf(msg, arg))
        else
            super.printf(msg, arg)
    }

    @Override
    void printf(String msg, Object[] args) {
        if( session?.quiet )
            return

        if( session?.ansiLog )
            log.info(String.printf(msg, args))
        else
            super.printf(msg, args)
    }

}
