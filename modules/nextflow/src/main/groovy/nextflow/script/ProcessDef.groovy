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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Const
import nextflow.Global
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.processor.TaskProcessor
import nextflow.script.dsl.ProcessConfigBuilder
import nextflow.script.params.BaseInParam
import nextflow.script.params.BaseOutParam
import nextflow.script.params.EachInParam

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

    ProcessDef(BaseScript owner, String name, ProcessConfig config, BodyDef taskBody) {
        this.owner = owner
        this.simpleName = name
        this.processName = name
        this.baseName = name
        this.processConfig = config
        this.taskBody = taskBody
    }

    /**
     * The process implementation, i.e. its {@code script}/{@code exec} body. Its
     * {@link BodyDef#source} is the process identity {@code TaskHasher} folds into every task
     * hash, so it is also what identifies this process when it is referenced from elsewhere
     * (see {@code AgentDef.toolsFingerprint}).
     */
    BodyDef getTaskBody() { taskBody }

    static String stripScope(String str) {
        str.split(Const.SCOPE_SEP).last()
    }

    @Override
    ProcessDef clone() {
        def result = (ProcessDef)super.clone()
        result.@processConfig = processConfig.clone()
        result.@taskBody = taskBody?.clone()
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
        return "Process `$name` declares ${expected} ${expected == 1 ? 'input' : 'inputs'} but was called with ${actual} ${actual == 1 ? 'argument' : 'arguments'}"
    }

    @Override
    Object run(Object[] args) {
        // invoke process with legacy inputs/outputs
        if( processConfig instanceof ProcessConfigV1 )
            output = runV1(args, processConfig)

        // invoke process with typed inputs/outputs
        else if( processConfig instanceof ProcessConfigV2 )
            output = runV2(args, processConfig)

        // return process output
        return output
    }

    private ChannelOut runV1(Object[] args, ProcessConfigV1 config) {
        // get params
        final params = ChannelOut.spread(args)
        final declaredInputs = config.getInputs()
        final declaredOutputs = config.getOutputs()

        // sanity check
        if( params.size() != declaredInputs.size() )
            throw new ScriptRuntimeException(missMatchErrMessage(processName, declaredInputs.size(), params.size()))

        // set input channels
        for( int i=0; i<params.size(); i++ ) {
            final inParam = (declaredInputs[i] as BaseInParam)
            inParam.setFrom(params[i])
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
                final param = (declaredOutputs[i] as BaseOutParam)
                final topicName = param.channelTopicName
                if( topicName && feedbackChannels )
                    throw new IllegalArgumentException("Output topic conflicts with recursion feature - process `$processName` should not declare any output topic" )
                final ch = feedbackChannels
                        ? feedbackChannels[i]
                        : ( topicName ? CH.createTopicSource(topicName) : CH.create(singleton) )
                param.setInto(ch)
            }
        }

        // make a copy of the output list because execution can change it
        final output = new ChannelOut(declaredOutputs.clone())

        // start processor
        createTaskProcessor().run()

        // the result channels
        assert declaredOutputs.size()>0, "Process output should contains at least one channel"
        return output
    }

    private ChannelOut runV2(Object[] args0, ProcessConfigV2 config) {
        final args = ChannelOut.spread(args0)
        final declaredInputs = config.getInputs()
        final declaredOutputs = config.getOutputs()

        // validate arguments
        if( args.size() != declaredInputs.size() )
            throw new ScriptRuntimeException(missMatchErrMessage(processName, declaredInputs.size(), args.size()))

        // set inputs
        for( int i = 0; i < declaredInputs.size(); i++ )
            declaredInputs[i].setChannel(createSourceChannel(args[i]))

        // set outputs
        final singleton = declaredInputs.isSingleton()

        final feedbackChannels = getFeedbackChannels()
        if( feedbackChannels && feedbackChannels.size() != declaredOutputs.size() )
            throw new ScriptRuntimeException("Process `$processName` inputs and outputs do not have the same cardinality - Feedback loop is not supported"  )

        final channels = new LinkedHashMap<String,DataflowWriteChannel>()
        for( int i = 0; i < declaredOutputs.size(); i++ ) {
            final param = declaredOutputs[i]
            final ch = feedbackChannels ? feedbackChannels[i] : CH.create(singleton)
            param.setChannel(ch)
            channels.put(param.getName(), ch)
        }

        for( final topic : declaredOutputs.getTopics() ) {
            final ch = CH.createTopicSource(topic.getTarget())
            topic.setChannel(ch)
        }

        // start processor
        createTaskProcessor().run()

        return new ChannelOut(channels)
    }

    private DataflowReadChannel createSourceChannel(Object value) {
        if( value instanceof DataflowReadChannel || value instanceof DataflowBroadcast )
            return CH.getReadChannel(value)

        final result = CH.value()
        result.bind(value)
        return result
    }

    TaskProcessor createTaskProcessor() {
        return createTaskProcessor(session, owner, processName, simpleName, baseName, processConfig, taskBody)
    }

    protected void applyConfig() {
        final configProcessScope = (Map)session.config.process
        new ProcessConfigBuilder(processConfig).applyConfig(configProcessScope, baseName, simpleName, processName)
    }

    /**
     * Build a {@link TaskProcessor} from a process config + body. Extracted as a
     * static helper so non-{@code ProcessDef} lowerings (e.g. {@code AgentDef} on
     * the tool-free task path) can drive the standard {@link TaskProcessor} pipeline
     * without duplicating the executor/factory wiring.
     *
     * NOTE: these static helpers overlap with {@link ProcessFactory}, whose only
     * remaining role is to be overridden by tests. Both should be consolidated into
     * a single static factory; deferred to keep the agent changeset small.
     */
    static TaskProcessor createTaskProcessor(Session session, BaseScript owner,
            String processName, String simpleName, String baseName,
            ProcessConfig config, BodyDef body) {
        // apply process directives from config settings
        new ProcessConfigBuilder(config).applyConfig((Map)session.config.process, baseName, simpleName, processName)

        return createTaskProcessorResolved(session, owner, processName, config, body)
    }

    /**
     * Build a task processor from an already resolved config. This allows non-process
     * components such as agents to use the canonical task/executor pipeline without
     * inheriting the {@code process} configuration scope.
     */
    static TaskProcessor createTaskProcessorResolved(Session session, BaseScript owner,
            String processName, ProcessConfig config, BodyDef body) {

        // create executor for process
        final executor = session
            .executorFactory
            .getExecutor(processName, config, body, session)

        // create task processor for process
        return session
            .newProcessFactory(owner)
            .newTaskProcessor(processName, executor, config, body)
    }

}
