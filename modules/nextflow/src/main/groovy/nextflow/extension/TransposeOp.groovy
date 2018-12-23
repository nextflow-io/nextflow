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
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.util.CheckHelper
/**
 * Implements channel transpose operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TransposeOp {

    static private Map TRANSPOSE_PARAMS = [ by: [Integer,List], remainder: Boolean ]

    private DataflowReadChannel source

    private DataflowWriteChannel target

    private List<Integer> cols

    private boolean remainder

    TransposeOp(DataflowReadChannel source, Map params=null) {
        CheckHelper.checkParams('transpose', params, TRANSPOSE_PARAMS)
        this.source = source
        this.target = new DataflowQueue()
        this.cols = parseByParam(params?.by)
        this.remainder = params?.remainder as Boolean
    }

    /**
     * Only for test -- do not use
     */
    @PackageScope
    TransposeOp() {}

    private List<Integer> parseByParam(value) {
        if( value == null )
            return Collections.emptyList()
        if( value instanceof List )
            return value as List<Integer>
        else
            return [value as int]
    }

    DataflowQueue apply() {
        DataflowHelper.subscribeImpl(source, [onNext: this.&transpose, onComplete: this.&done])
        return (DataflowQueue)target
    }

    protected void transpose(item) {
        item instanceof List ? transpose0(item) : target.bind(item)
    }

    protected List<Integer> findIndexes(List tuple) {
        def c = tuple.size()
        def result = new ArrayList(c)
        for( int i=0; i<c; i++ ) {
            if( tuple[i] instanceof List )
                result << i
        }
        return result
    }

    protected List<Integer> checkIndexes(List tuple, List<Integer> indexes) {
        for( int i=0; i<indexes.size(); i++ ) {
            def p = indexes[i]
            def v = tuple[p]
            if( !(v instanceof List) )
                throw new IllegalArgumentException("Not a valid transpose element at index: $p -- Offending tuple: (${tuple.join(', ')})")
        }
        return indexes
    }

    /**
     * Implements the transpose logic
     *
     * @param tuple The current channel tuple
     */
    protected void transpose0(List tuple) {
        // find the indexes of list elements in the provided tuple emitted by a channel
        List<Integer> indexes = cols ? checkIndexes(tuple,cols) : findIndexes(tuple)

        if( !indexes ) {
            target.bind(tuple)
            return
        }

        // find the max size of list elements to transpose
        int max=-1
        for( int i=0; i<indexes.size(); i++ ) {
            def c = ((List)tuple[ indexes[i] ]).size()
            if( c>max ) max=c
        }

        for( int i=0; i<max; i++ ) {
            final c = tuple.size()
            final buffer = new ArrayList(c)
            for( int j=0; j<c; j++ ) {
                // `j` is used to iterate over the `tuple`
                if( j in indexes ) {
                    def list = (List)tuple[j]
                    if( i<list.size() ) {
                        buffer[j] = list[i]
                    }
                    else if( remainder ) {
                        buffer[j] = null
                    }
                    else
                        return
                }
                else {
                    buffer[j] = tuple[j]
                }
            }

            target.bind(buffer)
        }
    }

    protected void done(dummy) {
        target.bind(Channel.STOP)
    }
}
