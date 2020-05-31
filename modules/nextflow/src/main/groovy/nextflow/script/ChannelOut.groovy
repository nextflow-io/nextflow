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
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.DuplicateChannelNameException
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList
import static nextflow.ast.NextflowDSLImpl.OUT_PREFIX
/**
 * Models the output of a process or a workflow component returning
 * more than one output channels
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelOut implements List<DataflowWriteChannel> {

    private @Delegate List<DataflowWriteChannel> target

    private Map<String,DataflowWriteChannel> channels

    ChannelOut() {
        this.target = Collections.<DataflowWriteChannel>unmodifiableList(Collections.<DataflowWriteChannel>emptyList())
        this.channels = Collections.<String,DataflowWriteChannel>emptyMap()
    }

    ChannelOut(List<DataflowWriteChannel> c) {
        this.target = Collections.<DataflowWriteChannel>unmodifiableList(c)
        this.channels = Collections.<String,DataflowWriteChannel>emptyMap()
    }

    ChannelOut(Map<String,DataflowWriteChannel> channels) {
        this.target = Collections.unmodifiableList(new ArrayList<DataflowWriteChannel>(channels.values()))
        this.channels = Collections.unmodifiableMap(new LinkedHashMap<String,DataflowWriteChannel>(channels))
    }

    ChannelOut(OutputsList outs) {
        channels = new HashMap<>(outs.size())
        final onlyWithName = new ArrayList<DataflowWriteChannel>(outs.size())
        for( OutParam param : outs ) {
            final ch = param.getOutChannel()
            final name = param.channelEmitName
            onlyWithName.add(ch)
            if(name) {
                if(channels.containsKey(name)) throw new DuplicateChannelNameException("Output channel name `$name` is used more than one time")
                channels.put(name, ch)
            }
        }
        target = Collections.unmodifiableList(onlyWithName)
    }

    Set<String> getNames() { channels.keySet().findAll { !it.startsWith(OUT_PREFIX) }  }

    @Override
    def getProperty(String name) {
        if( this.channels.containsKey(name) ) {
            return channels.get(name)
        }
        else
            metaClass.getProperty(this,name)
    }

    /**
     * Helper method that `spread`
     *
     * @param args
     * @return
     */
    static List spread(Object[] args) {
        final result = new ArrayList(args.size()*2)
        for( int i=0; i<args.size(); i++ ) {
            if( args[i] instanceof ChannelOut ) {
                final list = (List)args[i]
                for( def el : list ) {
                    result.add(el)
                }
            }
            else {
                result.add(args[i])
            }
        }
        return result
    }

    static List spread(List args) {
        spread(args as Object[])
    }

    static Object[] spreadToArray(Object[] args) {
        if(!args)
            return args
        // check if needed
        boolean found=false
        for( short i=0; i<args.length && !found; i++)
            found = args[i] instanceof ChannelOut
        if(!found)
            return args
        spread(args).toArray()
    }
}
