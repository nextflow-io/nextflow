/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.processor.TaskProcessor
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script implements ExecutionContext {

    private static Object[] EMPTY_ARGS = [] as Object[]

    private Session session

    private ProcessFactory processFactory

    private TaskProcessor taskProcessor

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
     * Holds the configuration object which will used to execution the user tasks
     */
    protected Map getConfig() {
        log.warn "The access of `config` object is deprecated"
        session.getConfig()
    }

    /**
     * Access to the last *process* object -- only for testing purpose
     */
    @PackageScope
    TaskProcessor getTaskProcessor() { taskProcessor }

    /**
     * Enable disable task 'echo' configuration property
     * @param value
     */
    protected void echo(boolean value = true) {
        log.warn "The use of `echo` method is deprecated"
        session.getConfig().process.echo = value
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
        binding.setVariable('moduleDir', meta.scriptPath?.parent )
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
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.preview.dsl=2` to allow the definition of workflow components")

        // launch the execution
        final workflow = new WorkflowDef(this, workflowBody)
        if( !binding.entryName )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected workflow(String name, Closure<BodyDef> workflowDef) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.preview.dsl=2` to allow the definition of workflow components")

        final workflow = new WorkflowDef(this,workflowDef,name)
        if( binding.entryName==name )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.preview.dsl=2` to import module files")

        include .setSession(session)
    }

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

        if( binding.entryName && !entryFlow ) {
            def msg = "Unknown workflow entry name: ${binding.entryName}"
            final allNames = meta.getLocalWorkflowNames()
            final guess = allNames.closest(binding.entryName)
            if( guess )
                msg += " -- Did you mean?\n" + guess.collect { "  $it"}.join('\n')
            throw new IllegalArgumentException(msg)
        }

        if( !entryFlow ) {
            log.debug "No entry workflow defined"
            return result
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        final ret = entryFlow.invoke_a(EMPTY_ARGS)
        session.notifyAfterWorkflowExecution()
        return ret
    }

    Object run() {
        setup()
        ExecutionStack.push(this)
        try {
            NF.dsl2 ? runDsl2() : runDsl1()
        }
        catch(InvocationTargetException e) {
            // provide the exception cause which is more informative than InvocationTargetException
            throw(e.cause ?: e)
        }
        finally {
            ExecutionStack.pop()
        }
    }

    protected abstract Object runScript()

}
