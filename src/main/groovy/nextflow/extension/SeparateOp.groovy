/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
