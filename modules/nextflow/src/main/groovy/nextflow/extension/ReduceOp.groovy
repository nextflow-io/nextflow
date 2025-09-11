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
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.ContextGrouping
import nextflow.extension.op.Op
/**
 * Implements reduce operator logic
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@CompileStatic
class ReduceOp {

    private DataflowReadChannel source
    private DataflowVariable target
    private Object seed
    private Closure action
    private Closure beforeBind

    private ReduceOp() { }

    ReduceOp(DataflowReadChannel source) {
        this.source = source
    }

    static ReduceOp create() {
        new ReduceOp()
    }

    ReduceOp withSource(DataflowReadChannel source) {
        this.source = source
        return this
    }

    ReduceOp withTarget(DataflowVariable target) {
        this.target = target
        return this
    }

    ReduceOp withSeed(Object seed) {
        this.seed = seed
        return this
    }

    ReduceOp withAction(Closure action) {
        this.action = action
        return this
    }

    ReduceOp withBeforeBind(Closure beforeBind) {
        this.beforeBind = beforeBind
        return this
    }

    private Session getSession() {
        return Global.session as Session
    }

    DataflowVariable apply() {
        if( source==null )
            throw new IllegalArgumentException("Missing reduce operator source channel")
        if( target==null )
            target = new DataflowVariable()
        final stopOnFirst = source instanceof DataflowExpression
        // the *accumulator* value
        def accum = this.seed

        // intercepts operator events
        final listener = new DataflowEventAdapter() {
            /*
             * when terminates bind the result value
             */
            void afterStop(final DataflowProcessor dp) {
                final result = beforeBind
                    ? beforeBind.call(accum)
                    : accum
                Op.bind(dp, target, result)
            }

            boolean onException(final DataflowProcessor dp, final Throwable e) {
                log.error("@unknown", e)
                session.abort(e)
                return true;
            }
        }

        final code = {
            final value = accum == null ? it : action.call(accum, it)
            final proc = getDelegate() as DataflowProcessor
            if( value!=Channel.VOID && value!=Channel.STOP ) {
                accum = value
            }
            if( stopOnFirst || value==Channel.STOP )
                proc.terminate()
        }

        new Op()
            .withInput(source)
            .withContext(new ContextGrouping())
            .withListener(listener)
            .withCode(code)
            .apply()

        return target
    }

}
