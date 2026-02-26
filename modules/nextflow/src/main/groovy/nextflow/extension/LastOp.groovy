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
 *
 */

package nextflow.extension


import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.extension.op.Op

/**
 * Implements last operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LastOp {

    private DataflowReadChannel source
    private DataflowVariable target

    LastOp withSource(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    LastOp withTarget(DataflowVariable target) {
        assert target!=null
        this.target = target
        return this
    }

    DataflowVariable apply() {
        assert source!=null
        if( target==null )
            target = new DataflowVariable()

        def last = null
        new SubscribeOp()
            .withInput(source)
            .withOnNext{ last = it }
            .withOnComplete{ DataflowProcessor dp -> Op.bind(dp, target, last) }
            .apply()
        return target
    }

}
