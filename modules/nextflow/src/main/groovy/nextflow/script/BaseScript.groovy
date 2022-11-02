/*
 * Copyright 2020-2022, Seqera Labs
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
import nextflow.exception.AbortOperationException
import nextflow.processor.TaskProcessor
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script implements ExecutionContext {

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
     * Access to the last *process* object -- only for testing purpose
     */
    @PackageScope
    TaskProcessor getTaskProcessor() { taskProcessor }

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

    protected process(String name, Closure<BodyDef> body) {
        def process = new ProcessDef(this,body,name)
        meta.addDefinition(process)
    }

    /**
     * Workflow main entry point
     *
     * @param body The implementation body of the workflow
     * @return The result of workflow execution
     */
    protected workflow(Closure<BodyDef> body) {
        final workflow = new WorkflowDef(this, body)
        if( !binding.entryName )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected workflow(String name, Closure<BodyDef> body) {
        final workflow = new WorkflowDef(this,body,name)
        if( binding.entryName==name )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(ExecutionStack.withinWorkflow())
            throw new IllegalStateException("Include statement is not allowed within a workflow definition")
        include .setSession(session)
    }

    @Override
    Object invokeMethod(String name, Object args) {
        binding.invokeMethod(name, args)
    }

    Object run() {
        setup()
        ExecutionStack.push(this)
        try {
            run0()
        }
        catch(InvocationTargetException e) {
            // provide the exception cause which is more informative than InvocationTargetException
            throw(e.cause ?: e)
        }
        finally {
            ExecutionStack.pop()
        }
    }

    private run0() {
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
            if( meta.getLocalWorkflowNames() )
                log.warn "No entry workflow specified"
            return result
        }

        // invoke the entry workflow
        session.notifyBeforeWorkflowExecution()
        final ret = entryFlow.invoke_a(BaseScriptConsts.EMPTY_ARGS)
        session.notifyAfterWorkflowExecution()
        return ret
    }

    protected abstract Object runScript()

    @Override
    void print(Object object) {
        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.print(object)
    }

    @Override
    void println() {
        if( session?.ansiLog )
            log.info("")
        else
            super.println()
    }

    @Override
    void println(Object object) {
        if( session?.ansiLog )
            log.info(object?.toString())
        else
            super.println(object)
    }

    @Override
    void printf(String msg, Object arg) {
        if( session?.ansiLog )
            log.info(String.printf(msg, arg))
        else
            super.printf(msg, arg)
    }

    @Override
    void printf(String msg, Object[] args) {
        if( session?.ansiLog )
            log.info(String.printf(msg, args))
        else
            super.printf(msg, args)
    }

}
