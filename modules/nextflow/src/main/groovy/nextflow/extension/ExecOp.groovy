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
import nextflow.script.ProcessDef
import nextflow.script.ChannelOut

/**
 * Implements the {@link OperatorImpl#exec} operator
 *
 * @author Jacob Munro <jacob.e.munro@gmail.com>
 */
@CompileStatic
class ExecOp {

    private ProcessDef processDef

    private DataflowReadChannel source

    private Object[] args 

    ExecOp( DataflowReadChannel source, ProcessDef processDef, Object[] args ) {
        assert processDef != null
        assert source != null

        this.source = source
        this.processDef = processDef
        this.args = args
    }

    Object apply() {

        def result = processDef.run(resolveInputs())

        return result
    }

    private Object[] resolveInputs() {

        Object[] resultArray = ([source] + (args as List)).toArray()

        return ChannelOut.spread(resultArray).toArray()
    }
}
