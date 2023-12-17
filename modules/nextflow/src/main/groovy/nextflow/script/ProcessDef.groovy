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
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Const
import nextflow.Global
import nextflow.NF
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.CombineManyOp
import nextflow.script.dsl.ProcessConfigBuilder

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
    private ProcessConfig config

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
        this.config = config
    }

    static String stripScope(String str) {
        str.split(Const.SCOPE_SEP).last()
    }

    protected void initialize() {
        // apply config settings to the process
        new ProcessConfigBuilder(config).applyConfig((Map)session.config.process, baseName, simpleName, processName)
    }

    @Override
    ProcessDef clone() {
        def result = (ProcessDef)super.clone()
        result.@taskBody = taskBody.clone()
        result.@config = config.clone()
        return result
    }

    @Override
    ProcessDef cloneWithName(String name) {
        ScriptMeta.addResolvedName(name)
        def result = clone()
        result.@processName = name
        result.@simpleName = stripScope(name)
        result.@config.processName = name
        return result
    }

    private ProcessInputs getDeclaredInputs() { config.getInputs() }

    private ProcessOutputs getDeclaredOutputs() { config.getOutputs() }

    BaseScript getOwner() { owner }

    String getName() { processName }

    String getSimpleName() { simpleName }

    String getBaseName() { baseName }

    ProcessConfig getProcessConfig() { config }

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

    private DataflowReadChannel collectInputs(Object[] args0) {
        final args = ChannelOut.spread(args0)
        if( args.size() != declaredInputs.size() )
            throw new ScriptRuntimeException(missMatchErrMessage(processName, declaredInputs.size(), args.size()))

        // emit value channel if process has no inputs
        if( args.size() == 0 ) {
            final source = CH.value()
            source.bind([])
            return source
        }

        // create input channels
        for( int i = 0; i < declaredInputs.size(); i++ )
            declaredInputs[i].bind(args[i])

        // combine input channels
        final inputs = declaredInputs.getChannels()
        if( inputs.size() == 1 )
            return inputs.first().chainWith( it -> [it] )

        final count = (0..<inputs.size()).count(
            i -> CH.isChannelQueue(inputs[i]) && !declaredInputs[i].isIterator()
        )
        if( NF.isStrictMode() && count > 1 )
            throw new ScriptRuntimeException("Process `$name` received multiple queue channel inputs which will be implicitly mergeed -- consider combining them explicitly with `combine` or `join`, or converting single-item chennels into value channels with `collect` or `first`")

        final iterators = (0..<inputs.size()).findAll( i -> declaredInputs[i].isIterator() )
        return CH.getReadChannel(new CombineManyOp(inputs, iterators).apply())
    }

    private void collectOutputs(boolean singleton) {
        // emit stdout if no outputs are defined
        if( declaredOutputs.size() == 0 ) {
            declaredOutputs.setDefault()
            return
        }

        // check for feedback channels
        final feedbackChannels = getFeedbackChannels()
        if( feedbackChannels && feedbackChannels.size() != declaredOutputs.size() )
            throw new ScriptRuntimeException("Process `$processName` inputs and outputs do not have the same cardinality - Feedback loop is not supported"  )

        for( int i=0; i<declaredOutputs.size(); i++ ) {
            final param = declaredOutputs[i]
            final topicName = param.getTopic()
            if( topicName && feedbackChannels )
                throw new IllegalArgumentException("Output topic conflicts with recursion feature - process `$processName` should not declare any output topic" )
            final ch = feedbackChannels
                    ? feedbackChannels[i]
                    : topicName ? CH.createTopicSource(topicName) : CH.create(singleton)
            param.setChannel(ch)
        }
    }

    @Override
    Object run(Object[] args) {
        // initialise process config
        initialize()

        // create input channel
        final source = collectInputs(args)

        // set output channels
        // note: the result object must be an array instead of a List to allow process
        // composition ie. to use the process output as the input in another process invocation
        final singleton = !CH.isChannelQueue(source)
        collectOutputs(singleton)

        // create the executor
        final executor = session
                .executorFactory
                .getExecutor(processName, config, taskBody, session)

        // create processor class
        session
                .newProcessFactory(owner)
                .newTaskProcessor(processName, executor, config, taskBody)
                .run(source)

        // the result channels
        return output = new ChannelOut(declaredOutputs)
    }

}
