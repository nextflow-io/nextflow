/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.script.ChannelOut

/**
 * Implements the {@link OperatorImpl#eval} operator
 *
 * @author Jacob Munro <jacob.e.munro@gmail.com>
 */
@CompileStatic
class EvalOp {

    private Object source
    
    private Closure closure

    EvalOp( ChannelOut source, Closure closure ) {
        assert source != null
        assert closure != null

        this.source = source
        this.closure = closure
    }

    EvalOp( DataflowReadChannel source, Closure closure ) {
        assert source != null
        assert closure != null

        this.source = source as DataflowWriteChannel
        this.closure = closure
    }

    Object apply() {

        final copy = (Closure)closure.clone()
        copy.setResolveStrategy(Closure.DELEGATE_FIRST)
        def result = source.with(copy)

        return result
    }

}
