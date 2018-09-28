/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.extension

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.SeparationClosure

/**
 * Implements the {@link DataflowExtensions#separate} operator logic
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SeparateOp {

    private DataflowReadChannel source

    private List<DataflowQueue> outputs

    private Closure mapper

    SeparateOp( final DataflowReadChannel source, final List<DataflowWriteChannel<?>> outputs, final Closure<List<Object>> code = null ) {
        assert source
        assert outputs

        this.source = source
        this.outputs = outputs.collect { (DataflowQueue)it }
        this.mapper = code

    }

    SeparateOp( final DataflowReadChannel source, final int n, Closure<List<Object>> mapper = null ) {
        assert source
        assert n

        this.source = source
        this.outputs = new ArrayList<>(n)
        for( int i=0; i<n; i++ ) {
            this.outputs[i] = new DataflowQueue()
        }

        this.mapper = mapper
    }

    @CompileDynamic
    private Closure<List<Object>> createDefaultMapper(int size) {

        int count = 0
        Closure<List<Object>> result = { it ->
            def tuple = it instanceof List ? it : [it]
            if( tuple.size() == size )
                return tuple

            else {
                if( count++ == 0 )
                    log.warn "The target channels number ($size) for the 'into' operator do not match the items number (${tuple.size()}) of the receveid tuple: $tuple"

                def result = new ArrayList(size)
                for( int i=0; i<size; i++ ) {
                    result[i] = i < tuple.size() ? tuple[i] : null
                }
                return result
            }
        }

        return result
    }

    List<DataflowQueue> apply() {
        if( !mapper )
            mapper = createDefaultMapper(outputs.size())

        DataflowHelper.newOperator( [source], outputs, new SeparationClosure(mapper))
        return outputs
    }


}
