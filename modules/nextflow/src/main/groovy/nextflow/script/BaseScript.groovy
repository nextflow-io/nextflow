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

import java.lang.reflect.InvocationTargetException
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.ast.ProcessFn
import nextflow.exception.AbortOperationException
import nextflow.script.dsl.ProcessDsl
import nextflow.script.dsl.WorkflowDsl
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseScript extends Script implements ExecutionContext {

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

    /**
     * Define a process.
     *
     * @param name
     * @param rawBody
     */
    protected void process(String name, Closure<BodyDef> rawBody) {
        final builder = new ProcessDsl(this, name)
        final copy = (Closure<BodyDef>)rawBody.clone()
        copy.delegate = builder
        copy.resolveStrategy = Closure.DELEGATE_FIRST
        final taskBody = copy.call()
        final process = builder.withBody(taskBody).build()
        meta.addDefinition(process)
    }

    /**
     * Define an anonymous workflow.
     *
     * @param rawBody
     */
    protected void workflow(Closure<BodyDef> rawBody) {
        final builder = new WorkflowDsl(this)
        final copy = (Closure<BodyDef>)rawBody.clone()
        copy.delegate = builder
        copy.resolveStrategy = Closure.DELEGATE_FIRST
        final body = copy.call()
        final workflow = builder.withBody(body).build()
        this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    /**
     * Define a named workflow.
     *
     * @param name
     * @param rawBody
     */
    protected void workflow(String name, Closure<BodyDef> rawBody) {
        final builder = new WorkflowDsl(this, name)
        final copy = (Closure<BodyDef>)rawBody.clone()
        copy.delegate = builder
        copy.resolveStrategy = Closure.DELEGATE_FIRST
        final body = copy.call()
        final workflow = builder.withBody(body).build()
        meta.addDefinition(workflow)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(ExecutionStack.withinWorkflow())
            throw new IllegalStateException("Include statement is not allowed within a workflow definition")
        include .setSession(session)
    }

    @CompileStatic
    private void registerProcessFunctions() {
        final dslEval = { Object delegate, Class<Closure> clazz ->
            final cl = clazz.newInstance(this, this)
            cl.delegate = delegate
            cl.resolveStrategy = Closure.DELEGATE_FIRST
            cl.call()
        }

        final clazz = this.getClass()
        for( final method : clazz.getDeclaredMethods() ) {
            // check for ProcessFn annotation
            final name = method.getName()
            final processFn = method.getAnnotation(ProcessFn)
            if( !processFn )
                continue

            // validate annotation
            if( processFn.script() && processFn.shell() )
                throw new IllegalArgumentException("Process function `${name}` cannot have script and shell enabled simultaneously")

            // build process from annotation
            final builder = new ProcessDsl(this, name)
            dslEval(builder, processFn.directives())
            dslEval(builder, processFn.inputs())
            dslEval(builder, processFn.outputs())

            // get method parameters
            final paramNames = (List<String>)((Closure)processFn.params().newInstance(this, this)).call()
            final params = (0 ..< paramNames.size()).collect( i ->
                new Parameter( paramNames[i], method.getParameters()[i].getType() )
            )

            builder.config.params = params

            // determine process type
            def type
            if( processFn.script() )
                type = 'script'
            else if( processFn.shell() )
                type = 'shell'
            else
                type = 'exec'

            // create task body
            final taskBody = new BodyDef( this.&"${name}", processFn.source(), type, [] )

            final process = builder.withBody(taskBody).build()
            meta.addDefinition(process)
        }
    }

    private run0() {
        // register any process functions
        registerProcessFunctions()

        // execute script
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
