/*
 * Copyright 2020-2021, Seqera Labs
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

import nextflow.script.testflow.TestCase
import nextflow.script.testflow.TestFailure
import nextflow.script.testflow.TestSuite

import java.lang.reflect.InvocationTargetException
import java.nio.file.Paths

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.processor.TaskProcessor

import java.time.Duration
import java.time.Instant
import java.time.LocalDateTime
import java.time.temporal.ChronoUnit

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

    private List<TestflowDef> testFlows = new ArrayList<>()

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
        final msg = "The access of `config` object is deprecated"
        if( NF.dsl2Final )
            throw new DeprecationException(msg)
        log.warn(msg)
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
        final msg = "The use of `echo` method has been deprecated"
        if( NF.dsl2Final )
            throw new DeprecationException(msg)
        log.warn(msg)
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
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.enable.dsl=2` to allow the definition of workflow components")

        // launch the execution
        final workflow = new WorkflowDef(this, workflowBody)
        if( !binding.entryName )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected workflow(String name, Closure<BodyDef> workflowDef) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.enable.dsl=2` to allow the definition of workflow components")

        final workflow = new WorkflowDef(this,workflowDef,name)
        if( binding.entryName==name )
            this.entryFlow = workflow
        meta.addDefinition(workflow)
    }

    protected IncludeDef include( IncludeDef include ) {
        if(!NF.isDsl2())
            throw new IllegalStateException("Module feature not enabled -- Set `nextflow.enable.dsl=2` to import module files")

        include .setSession(session)
    }

    protected testflow(String name, Closure testflowBody) {
        if( !name )
            throw new IllegalArgumentException("Missing testflow name")
        final test = new TestflowDef(this, name, testflowBody)
        testFlows << test
        meta.addDefinition(test)
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

        if( binding.testMode ) {
            if( !testFlows ) throw new IllegalArgumentException("No tests defined")
            for( TestflowDef test : testFlows ) {
                test.invoke_a(EMPTY_ARGS)
            }
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

    TestSuite checkTests() {
        if (!binding.testMode)
            throw new IllegalStateException("Not running in test mode")

        int tests = 0
        int skipped = 0
        int failures = 0
        int errors = 0
        final scriptName = session.testName ?: session.scriptName.replace(".nf", "")

        Duration suiteTime = Duration.ZERO
        List<TestCase> testCase = new ArrayList<>()
        for (TestflowDef test : testFlows) {

            tests += 1
            suiteTime += test.runTime

            TestCase testcase = new TestCase(name: test.name, className: scriptName, time: test.getRunTime())
            try {
                test.validateExecution()
            } catch (AssertionError e) {
                failures += 1
                testcase.failure = new TestFailure(
                        message: "assertion error",
                        type: e.class.name,
                        content: e.message
                )

                if (e instanceof TestflowDsl.TaskAssertionError) {
                    testcase.workDir = e.context.meta.workDir
                }

            } catch (RuntimeException e) {
                errors += 1
                testcase.failure = new TestFailure(
                        message: "runtime error",
                        type: e.class.name,
                        content: e.message
                )
            }

            testCase.add(testcase)

        }

        return new TestSuite(
                name: scriptName,
                tests: tests,
                skipped: skipped,
                failures: failures,
                errors: errors,
                testcase: testCase,
                systemErr: "",
                systemOut: "",
                time: suiteTime,
                timestamp: LocalDateTime.now().truncatedTo(ChronoUnit.SECONDS)
        )
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
