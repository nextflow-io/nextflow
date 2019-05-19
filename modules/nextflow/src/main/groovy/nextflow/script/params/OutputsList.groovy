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

import groovyx.gpars.dataflow.DataflowWriteChannel


/**
 * Container to hold all process outputs
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OutputsList implements List<OutParam>, Cloneable {

    @Override
    OutputsList clone() {
        def result = (OutputsList)super.clone()
        result.target = new ArrayList<>(target.size())
        for( OutParam param : target )
            result.add((OutParam)param.clone())
        return result
    }

    @Delegate
    private List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() {
        final List<DataflowWriteChannel> result = new ArrayList<>(target.size())
        for(OutParam param : target) { result.addAll(param.getOutChannels()) }
        return result
    }

    List<String> getNames() { target *. name }

    def <T extends OutParam> List<T> ofType( Class<T>... classes ) {
        (List<T>) target.findAll { it.class in classes }
    }

    void setSingleton( boolean value ) {
        for( OutParam param : target ) { ((BaseOutParam)param).singleton = value }
    }
}
