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

package nextflow.script.params

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.extension.CH


/**
 * Container to hold all process outputs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class InputsList implements List<InParam>, Cloneable {

    @Override
    InputsList clone() {
        def result = (InputsList)super.clone()
        result.target = new ArrayList<>(target.size())
        for( InParam param : target ) {
            result.target.add((InParam)param.clone())
        }
        return result
    }

    @Delegate
    private List<InParam> target = new LinkedList<>()

    List<DataflowReadChannel> getChannels() {
        target.collect { InParam it -> it.getInChannel() }
    }

    List<String> getNames() { target *. name }


    def <T extends InParam> List<T> ofType( Class<T> clazz ) {
        (List<T>) target.findAll { it.class == clazz }
    }

    boolean allScalarInputs() {
        for( InParam param : target ) {
            if( CH.isChannelQueue(param.inChannel) )
                return false
        }
        return true
    }

}

