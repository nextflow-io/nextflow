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
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.Op

/**
 * Implements "view" operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ViewOp {

    private DataflowReadChannel source
    private Closure<Boolean> action
    private boolean newLine

    private static Session getSession() { Global.getSession() as Session }

    ViewOp() {
        this.source = source
        this.action = action
    }

    ViewOp withSource(DataflowReadChannel source) {
        if( source==null )
            throw new IllegalArgumentException("Argument 'source' cannot be null")
        this.source = source
        return this
    }

    ViewOp withCode(Closure<Boolean> action) {
        this.action = action
        return this
    }

    ViewOp withNewLine(boolean newLine) {
       this.newLine = newLine
        return this
    }

    DataflowWriteChannel apply() {
        final target = CH.createBy(source);

        new SubscribeOp()
            .withInput(source)
            .withOnNext { DataflowProcessor dp, Object it->
                final obj = action != null ? action.call(it) : it
                session.printConsole(obj?.toString(), newLine)
                Op.bind(dp,target,it)
            }
            .withOnComplete { CH.close0(target) }
            .apply()

        return target
    }
}
