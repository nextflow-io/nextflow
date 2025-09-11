/*
 * Copyright 2013-2024, Seqera Labs
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


import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.extension.op.Op
import nextflow.splitter.FastqSplitter
/**
 * Implements the inner merging logic for splitterXxx operator(s)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class SplitterMergeClosure extends Closure {

    private int numOfParams

    private List<Integer> indexes

    private int emissionCount

    private DataflowWriteChannel target

    SplitterMergeClosure(List<Integer> indexes, DataflowWriteChannel target) {
        super(null, null);
        this.numOfParams = indexes.size()
        this.indexes = indexes
        this.target = target
    }

    @Override
    int getMaximumNumberOfParameters() {
        numOfParams
    }

    @Override
    Class[] getParameterTypes() {
        Collections.nCopies(numOfParams, Object)
    }

    @Override
    void setDelegate(final Object delegate) {
        super.setDelegate(delegate);
    }

    @Override
    void setResolveStrategy(final int resolveStrategy) {
        super.setResolveStrategy(resolveStrategy);
    }

    @Override
    Object call(final Object arguments) {
        throw new UnsupportedOperationException()
    }

    @Override
    Object call(final Object... args) {
        log.trace "merging ($emissionCount) indexes=$indexes; args=$args"
        List result = null
        boolean header = false
        for( int i=0; i<args.size(); i++ ) {
            final item = args[i]
            // - When Pair-ended splitting is enabled the Fastq splitter
            //   emits the indexes of the split files as first tuple
            // - Those indexes are needed to merge the final tuple in correct order
            if( emissionCount==0 && item instanceof FastqSplitter.SplitIndex ) {
                indexes[i] = item.value
                header |= true
            }

            // - Compose the i-th tuple
            else if( item instanceof List ) {
                List entry = item
                if( i==0 ) {
                    result = new ArrayList(entry)
                }
                result[indexes[i]] = entry[indexes[i]]
            }
            else
                throw new IllegalArgumentException("Invalid splitter entry -- offending value=$item; indexes=$indexes; args=$args")

        }

        // emit the merged tuple (skipping the header)
        final dp = getDelegate() as DataflowProcessor
        if( !header && dp ) {
            log.trace "merging ($emissionCount) result=$result"
            Op.bind(dp, target, result)
        }

        emissionCount++

        return result;
    }

    @Override
    Object call() {
        throw new UnsupportedOperationException()
    }
}
