/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.script.params.v2

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.script.params.OutParam

/**
 * Models the process outputs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessOutputs implements Cloneable {

    private List<ProcessOutput> params = []

    private List<ProcessTopic> topics = []

    /**
     * Environment variables which will be exported from the
     * task environment for each task and made available to
     * process outputs.
     */
    private Set<String> env = []

    /**
     * Shell commands which will be executed in the task environment
     * for each task and whose output will be made available
     * to process outputs. The key corresponds to the environment
     * variable to which the command output will be saved.
     */
    private Map<String,Object> eval = [:]

    /**
     * Output files which will be unstaged from the task
     * directory for each task and made available to process
     * outputs.
     */
    private Map<String,ProcessFileOutput> files = [:]

    void addParam(String name, Class type, Object value) {
        params.add(new ProcessOutput(name, type, value))
    }

    void addTopic(Object value, String target) {
        topics.add(new ProcessTopic(value, target))
    }

    void addEnv(String name) {
        env.add(name)
    }

    void addEval(String name, Object value) {
        eval.put(name, value)
    }

    void addFile(String key, ProcessFileOutput file) {
        files.put(key, file)
    }

    List<ProcessOutput> getParams() {
        return params
    }

    int size() {
        return params.size()
    }

    ProcessOutput getAt(int i) {
        return params.get(i)
    }

    List<ProcessTopic> getTopics() {
        return topics
    }

    Set<String> getEnv() {
        return env
    }

    Map<String,Object> getEval() {
        return eval
    }

    Map<String,ProcessFileOutput> getFiles() {
        return files
    }

    List<DataflowWriteChannel> getChannels() {
        final result = new HashSet<DataflowWriteChannel>()
        for( final param : params )
            result.add(param.getOutChannel())
        for( final topic : topics )
            result.add(topic.getChannel())
        return new ArrayList<>(result)
    }

    @Override
    ProcessOutputs clone() {
        final result = (ProcessOutputs)super.clone()
        result.params = new ArrayList<>(params.size())
        for( final param : params )
            result.params.add(param.clone())
        return result
    }

}


@CompileStatic
class ProcessOutput implements OutParam {

    /**
     * Name of the output channel in the process outputs (i.e. `.out`).
     */
    private String name

    /**
     * Optional output type which is used to validate task outputs.
     */
    private Class type

    /**
     * Lazy expression (e.g. closure) which defines the output value
     * in terms of the task context, including environment variables,
     * files, and standard output.
     * It will be evaluated for each task after it is executed. 
     */
    private Object value

    /**
     * Output channel which is created when the process is invoked
     * in a workflow.
     */
    DataflowWriteChannel channel

    ProcessOutput(String name, Class type, Object value) {
        this.name = name
        this.type = type
        this.value = value
    }

    @Override
    String getName() {
        return name
    }

    Object getLazyValue() {
        return value
    }

    @Override
    ProcessOutput clone() {
        (ProcessOutput)super.clone()
    }

    /// LEGACY METHODS

    DataflowWriteChannel getOutChannel() {
        return channel
    }

    short getIndex() {
        throw new UnsupportedOperationException()
    }

    String getChannelEmitName() {
        throw new UnsupportedOperationException()
    }

    String getChannelTopicName() {
        throw new UnsupportedOperationException()
    }

}


@CompileStatic
class ProcessTopic {

    /**
     * Lazy expression (e.g. closure) which defines the output value
     * in terms of the task context, including environment variables,
     * files, and standard output.
     * It will be evaluated for each task after it is executed. 
     */
    private Object value

    /**
     * Name of the target topic.
     */
    private String target

    /**
     * Topic channel which is bound when the process is invoked
     * in a workflow.
     */
    DataflowWriteChannel channel

    ProcessTopic(Object value, String target) {
        this.value = value
        this.target = target
    }

    Object getLazyValue() {
        return value
    }

    String getTarget() {
        return target
    }

    @Override
    ProcessTopic clone() {
        (ProcessTopic)super.clone()
    }

}
