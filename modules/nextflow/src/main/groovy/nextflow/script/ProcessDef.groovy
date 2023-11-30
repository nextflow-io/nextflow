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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Const
import nextflow.Global
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.MergeOp
import nextflow.script.dsl.ProcessBuilder
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
class ProcessDef extends BindableDef implements IterableDef, ChainableDef {

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
     * The resolved process configuration
     */
    private ProcessConfig processConfig

    /**
     * The actual process implementation
     */
    private BodyDef taskBody

    /**
     * The result of the process execution
     */
    private transient ChannelOut output

    ProcessDef(BaseScript owner, String name, BodyDef body, ProcessConfig config) {
        this.owner = owner
        this.simpleName = name
        this.processName = name
        this.baseName = name
        this.taskBody = body
        this.processConfig = config
    }

    static String stripScope(String str) {
        str.split(Const.SCOPE_SEP).last()
    }

    protected void initialize() {
        // apply config settings to the process
        new ProcessBuilder(processConfig).applyConfig((Map)session.config.process, baseName, simpleName, processName)
    }

    @Override
    ProcessDef clone() {
        def result = (ProcessDef)super.clone()
        result.@taskBody = taskBody.clone()
        result.@processConfig = processConfig.clone()
        return result
    }

    @Override
    ProcessDef cloneWithName(String name) {
        ScriptMeta.addResolvedName(name)
        def result = clone()
        result.@processName = name
        result.@simpleName = stripScope(name)
        result.@processConfig.processName = name
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
        if( output==null )
            throw new ScriptRuntimeException("Access to '${processName}.out' is undefined since the process '$processName' has not been invoked before accessing the output attribute")
        if( output.size()==0 )
            throw new ScriptRuntimeException("Access to '${processName}.out' is undefined since the process '$processName' doesn't declare any output")
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

        // create merged input channel
        final inputs = ChannelOut.spread(args).collect(ch -> getInChannel(ch))
        if( inputs.findAll(ch -> CH.isChannelQueue(ch)).size() > 1 )
            throw new ScriptRuntimeException("Process `$name` received multiple queue channel inputs which is not allowed")

        final input = CH.getReadChannel(new MergeOp(inputs.first(), inputs[1..<inputs.size()], [flat: false]).apply())

        // set input channels
        for( int i=0; i<declaredInputs.size(); i++ ) {
            final inParam = (declaredInputs[i] as BaseInParam)
            inParam.setFrom("in${i}")
            inParam.init()
        }

        // set output channels
        // note: the result object must be an array instead of a List to allow process
        // composition ie. to use the process output as the input in another process invocation
        if( declaredOutputs.size() ) {
            final allScalarValues = declaredInputs.allScalarInputs()
            final hasEachParams = declaredInputs.any { it instanceof EachInParam }
            final singleton = allScalarValues && !hasEachParams

            // check for feedback channels
            final feedbackChannels = getFeedbackChannels()
            if( feedbackChannels && feedbackChannels.size() != declaredOutputs.size() )
                throw new ScriptRuntimeException("Process `$processName` inputs and outputs do not have the same cardinality - Feedback loop is not supported"  )

            for(int i=0; i<declaredOutputs.size(); i++ ) {
                final ch = feedbackChannels ? feedbackChannels[i] : CH.create(singleton)
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
                .run(input)

        // the result channels
        assert declaredOutputs.size()>0, "Process output should contains at least one channel"
        return output = new ChannelOut(copyOuts)
    }

    private DataflowReadChannel getInChannel(Object obj) {
        if( obj == null )
            throw new IllegalArgumentException('A process input channel evaluates to null')

        final result = obj instanceof Closure
            ? obj.call()
            : obj

        if( result == null )
            throw new IllegalArgumentException('A process input channel evaluates to null')

        def inChannel
        if ( result instanceof DataflowReadChannel || result instanceof DataflowBroadcast )
            inChannel = CH.getReadChannel(result)
        else {
            inChannel = CH.value()
            inChannel.bind(result)
        }

        return inChannel
    }

}
