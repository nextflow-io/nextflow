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
import static nextflow.ast.NextflowDSLImpl.OUT_PREFIX
/**
 * Models the output of a process or a workflow component returning
 * more than one output channels
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelOut implements List {

    final private static Map<String,Integer> NUMS = new HashMap<>(10)

    static {
        NUMS['first'] = 0
        NUMS['second'] = 1
        NUMS['third'] = 2
        NUMS['fourth'] = 3
        NUMS['fifth'] = 4
        NUMS['sixth'] = 5
        NUMS['seventh'] = 6
        NUMS['eighth'] = 7
        NUMS['ninth'] = 8
        NUMS['tenth'] = 9
    }

    @Delegate List target

    private Map<String,Object> channels

    ChannelOut() {
        this.target = Collections.<Object>unmodifiableList(Collections.<Object>emptyList())
        this.channels = Collections.emptyMap()
    }

    ChannelOut(List c) {
        target = Collections.unmodifiableList(c)
        this.channels = Collections.emptyMap()
    }

    ChannelOut(LinkedHashMap<String,?> channels) {
        this.target = Collections.unmodifiableList(new ArrayList(channels.values()))
        this.channels = new LinkedHashMap(channels)
    }

    Set<String> getNames() { channels.keySet().findAll { !it.startsWith(OUT_PREFIX) }  }

    def getProperty(String name) {
        if( channels.containsKey(name) ) {
            return channels.get(name)
        }
        else if(NUMS.containsKey(name)) {
            // this has been deprecated
            final i = NUMS[name]
            log.warn1 "Property `out.${name}` has been deprecated -- Use `out[$i]` instead"
            return target[i]
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
}
