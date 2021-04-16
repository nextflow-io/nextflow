/*
 * Copyright 2020, Seqera Labs
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

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Channel
import nextflow.extension.CH
import nextflow.processor.TaskMeta

/**
 * Wrap a process or workflow output
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TestflowDsl {

    private int current
    private String type
    private String name
    private int outputsCount
    private List<ChannelEntry> channelEntries
    private DataflowQueue tasksMetaChannel
    @Lazy private List<EmissionValues> emissions = { fetchOutputs() }()
    @Lazy private List<TaskMeta> tasksMetaList = { fetchTasksMeta() }()

    TestflowDsl(ProcessDef process) {
        this.type = process.type
        this.name = process.name
        this.tasksMetaChannel = process.tasksMeta
        init(process.getOut())
    }

    TestflowDsl(WorkflowDef workflow) {
        this.type = workflow.type
        this.name = workflow.name
        init(workflow.getOut())
    }

    protected TestflowDsl(ChannelOut outputs) {
        init(outputs)
    }

    private void init(ChannelOut outputs) {
        this.outputsCount = outputs.size()
        this.channelEntries = readable0(outputs)
    }

    private List<ChannelEntry> readable0(ChannelOut output) {
        final result = new ArrayList(output.size())
        for( int i=0; i<output.size(); i++ ) {
            final ch = output.get(i)
            final name = output.nameOf(ch)
            final read = CH.getReadChannel(ch)
            result << new ChannelEntry(i, name, read)
        }
        return result
    }

    private List<TaskMeta> fetchTasksMeta() {
        if( tasksMetaChannel==null )
            return null
        final result = new ArrayList()
        while( true ) {
            final value = tasksMetaChannel.getVal()
            if( value == Channel.STOP )
                break
            result.add(value)
        }
        return result
    }

    private List<EmissionValues> fetchOutputs() {
        List<EmissionValues> result = []
        boolean single=true
        boolean terminated=false
        for( int i=0; !terminated && channelEntries.size()>0; i++ ) {
            def values = new ArrayList()
            def named = new HashMap()
            for( int j=0; !terminated && j<channelEntries.size(); j++ ) {
                final x = channelEntries[j].read()
                log.trace "Read test output value: '$x' produced by $type '$name'"
                if( Channel.STOP.is(x) ) {
                    terminated=true
                    break
                }
                values[j] = x
                if( channelEntries[j].name )
                    named[ channelEntries[j].name ] = x
                if( i==0 )
                    single &= channelEntries[j].isValue()
            }
            //
            if( !terminated ) {
                result << new EmissionValues(values, named, this)
            }
            if( single ) {
                // if channel are all dataflow variables stop at the first iteration
                break
            }
        }
        return result
    }

    int outputsCount() {
        return outputsCount
    }

    int emissionsCount() {
        return emissions.size()
    }

    void emissionWith(Map opts, Closure body) {
        if( type!='process' )
            throw new IllegalArgumentException("Test 'emissionWith' can only be applied to process outputs")

        int pos
        if( opts.tag!=null ) {
            pos = tasksMetaList.findIndexOf { it.tag == opts.tag }
            if( pos==-1 ) throw new IllegalArgumentException("Unable to find any task with tag=$opts.tag")
        }
        else if( opts.index!=null ) {
            pos = tasksMetaList.findIndexOf { it.index == opts.index }
            if( pos==-1 ) throw new IllegalArgumentException("Unable to find any task with index=$opts.index")
        }
        else
            throw new IllegalArgumentException("Test 'emissionWith' requires 'tag' or 'index' argument")

        apply0(pos,body)
    }

    void emissionNext(Closure body) {
        apply0(current++, body)
    }

    private void apply0(int index, Closure body) {
        if( type=='process' ) {
            final meta = tasksMetaList.get(index)
            final values = emissions.get(index)
            final ctx = new TaskDsl( values, meta )
            try {
                invoke(body, ctx)
            } catch(AssertionError e) {
                throw new TaskAssertionError(message: e.message,  context: ctx)
            }
        }
        else {
            final values = emissions.get(index)
            final ctx = new WorkflowDsl( values )
            try {
                invoke(body, ctx)
            } catch(AssertionError e) {
                throw new WorkflowAssertionError(message: e.message, context: ctx)
            }
        }
    }

    private void invoke(Closure closure, Object binding) {
        def copy = (Closure) closure.clone()
        copy.setDelegate(binding)
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        copy.call()
    }

    String toString() { return "$type:$name" }


    @TupleConstructor(includeFields = true)
    @ToString(includeNames = true, includePackage = false)
    static class ChannelEntry {
        private int index
        private String name
        private DataflowReadChannel channel

        int getIndex() { index }
        String getName() { name }
        boolean isValue() { channel instanceof DataflowExpression }
        def read() { channel.getVal() }
    }

    @TupleConstructor(includeFields = true)
    static class EmissionValues {
        private List values
        private Map<String,Object> named
        private TestflowDsl parent

        int size() { return values.size() }
        int getLength() { values.size() }
        Object getAt(int index) { values[index] }

        @Override
        def getProperty(String name) {
            this.named.containsKey(name)
                ? named.get(name)
                : metaClass.getProperty(this,name)
        }

        String toString() { return "${parent.type}:${parent.name}" }
    }

    @Canonical
    @ToString(includeNames = true, includePackage = false)
    static class TaskDsl {
        EmissionValues out
        TaskMeta meta

        Path path(String path) {
            meta.workDir.resolve(path)
        }
    }

    @Canonical
    static class TaskAssertionError extends AssertionError {
        String message
        TaskDsl context
    }

    @Canonical
    @ToString(includeNames = true, includePackage = false)
    static class WorkflowDsl {
        EmissionValues out
    }

    @Canonical
    static class WorkflowAssertionError extends AssertionError {
        String message
        WorkflowDsl context
    }


}
