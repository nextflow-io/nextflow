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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.Global
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.script.params.BaseInParam
import nextflow.script.params.BaseOutParam
import nextflow.script.params.EachInParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutputsList

/**
 * Models a nextflow process definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessDef extends BindableDef implements ChainableDef {

    private Session session = Global.session as Session

    /**
     * The script owning this process
     */
    private BaseScript owner

    /**
     * Fully qualified process name ie. it may contain nested workflow execution scopes
     */
    private String processName

    /**
     * Simple process names ie. the name used to declared or import it w/o the execution scope
     * This name is used to resolve the configuration
     */
    private String simpleName

    /**
     * Process source name how it was given in the original model ie. this name
     * name is stable and cannot be changed/aliased while the process is included.
     */
    private String baseName

    /**
     * The closure holding the process definition body
     */
    private Closure<BodyDef> rawBody

    /**
     * The resolved process configuration
     */
    private transient ProcessConfig processConfig

    /**
     * The actual process implementation
     */
    private transient BodyDef taskBody

    /**
     * The result of the process execution
     */
    private transient ChannelOut output

    ProcessDef(BaseScript owner, Closure<BodyDef> body, String name ) {
        this.owner = owner
        this.rawBody = body
        this.simpleName = name
        this.processName = name
        this.baseName = name
    }

    static String stripScope(String str) {
        str.split(Const.SCOPE_SEP).last()
    }

    protected void initialize() {
        log.trace "Process config > $processName"
        assert processConfig==null

        // the config object
        processConfig = new ProcessConfig(owner,processName)

        // Invoke the code block which will return the script closure to the executed.
        // As side effect will set all the property declarations in the 'taskConfig' object.
        processConfig.throwExceptionOnMissingProperty(true)
        final copy = (Closure)rawBody.clone()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.setDelegate(processConfig)
        taskBody = copy.call() as BodyDef
        processConfig.throwExceptionOnMissingProperty(false)
        if ( !taskBody )
            throw new ScriptRuntimeException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")

        // apply config settings to the process
        processConfig.applyConfig((Map)session.config.process, baseName, simpleName, processName)
    }

    @Override
    ProcessDef clone() {
        def result = (ProcessDef)super.clone()
        result.@taskBody = taskBody?.clone()
        result.@rawBody = (Closure)rawBody?.clone()
        return result
    }

    @Override
    ProcessDef cloneWithName(String name) {
        def result = clone()
        result.@processName = name
        result.@simpleName = stripScope(name)
        return result
    }

    private InputsList getDeclaredInputs() { processConfig.getInputs() }

    private OutputsList getDeclaredOutputs() { processConfig.getOutputs() }

    BaseScript getOwner() { owner }

    String getName() { processName }

    String getSimpleName() { simpleName }

    String getBaseName() { baseName }

    ProcessConfig getProcessConfig() { processConfig }

    ChannelOut getOut() {
        if(!output) throw new ScriptRuntimeException("Access to '${processName}.out' is undefined since process doesn't declare any output")
        return output
    }

    String getType() { 'process' }

    private String missMatchErrMessage(String name, int expected, int actual) {
        final ch = expected > 1 ? "channels" : "channel"
        return "Process `$name` declares ${expected} input ${ch} but ${actual} were specified"
    }

    @Override
    Object run(Object[] args) {
        // initialise process config
        initialize()

        // get params 
        final params = ChannelOut.spread(args)
        // sanity check
        if( params.size() != declaredInputs.size() )
            throw new ScriptRuntimeException(missMatchErrMessage(processName, declaredInputs.size(), params.size()))

        // set input channels
        for( int i=0; i<params.size(); i++ ) {
            (declaredInputs[i] as BaseInParam).setFrom(params[i])
        }

        // set output channels
        // note: the result object must be an array instead of a List to allow process
        // composition ie. to use the process output as the input in another process invocation
        if( declaredOutputs.size() ) {
            List result = new ArrayList<>(declaredOutputs.size())
            final allScalarValues = declaredInputs.allScalarInputs()
            final hasEachParams = declaredInputs.any { it instanceof EachInParam }
            final singleton = allScalarValues && !hasEachParams

            for(int i=0; i<declaredOutputs.size(); i++ ) {
                final ch = CH.create(singleton)
                result[i] = ch 
                (declaredOutputs[i] as BaseOutParam).setInto(ch)
            }
        }

        // make a copy of the output list because execution can change it
        final copyOuts = declaredOutputs.clone()

        // create the executor
        final executor = session
                .executorFactory
                .getExecutor(processName, processConfig, taskBody, session)

        // create processor class
        session
                .newProcessFactory(owner)
                .newTaskProcessor(processName, executor, processConfig, taskBody)
                .run()

        // the result channels
        assert declaredOutputs.size()>0, "Process output should contains at least one channel"
        return output = new ChannelOut(copyOuts)
    }

}
