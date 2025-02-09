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
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Global
import nextflow.Session
import nextflow.extension.op.Op
import nextflow.extension.op.OpContext
/**
 * Implements the "subscribe" operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SubscribeOp {

    static private List<String> VALID_HANDLERS = [ 'onNext', 'onComplete', 'onError' ]

    private DataflowReadChannel source
    private OpContext context
    private Closure onNext
    private Closure onComplete
    private Closure onError

    private static Session getSession() { Global.getSession() as Session }

    SubscribeOp withInput(DataflowReadChannel source) {
        assert source!=null
        this.source = source
        return this
    }

    SubscribeOp withOnNext(Closure event) {
        this.onNext = event
        return this
    }

    SubscribeOp withOnComplete(Closure event) {
        this.onComplete = event
        return this
    }

    SubscribeOp withOnError(Closure event) {
        this.onError = event
        return this
    }

    SubscribeOp withEvents(Map<String,Closure> events) {
        if( events ) {
            events.keySet().each {
                if( !VALID_HANDLERS.contains(it) )  throw new IllegalArgumentException("Not a valid handler name: $it")
            }
            this.onNext = events.onNext as Closure
            this.onComplete = events.onComplete as Closure
            this.onError = events.onError as Closure
        }
        return this
    }

    SubscribeOp withContext(OpContext context) {
        if( context!=null )
            this.context = context
        return this
    }

    DataflowProcessor apply() {
        def error = false
        def stopOnFirst = source instanceof DataflowExpression
        def listener = new DataflowEventAdapter() {

            @Override
            void afterStop(final DataflowProcessor dp) {
                if( !onComplete || error ) return
                try {
                    onComplete.call(dp)
                }
                catch( Exception e ) {
                    SubscribeOp.log.error("@unknown", e)
                    session.abort(e)
                }
            }

            @Override
            boolean onException(final DataflowProcessor dp, final Throwable e) {
                error = true
                if( !onError ) {
                    log.error("@unknown", e)
                    session.abort(e)
                }
                else {
                    onError.call(e)
                }
                return true
            }
        }

        final code = {
            final proc = getDelegate() as DataflowProcessor
            if( onNext instanceof Closure ) {
                final action = (Closure) onNext
                final types = action.getParameterTypes()
                types.size()==2 && types[0]==DataflowProcessor.class
                    ? action.call(proc, it)
                    : action.call(it)
            }
            if( stopOnFirst ) {
                proc.terminate()
            }
        }

        new Op()
            .withInput(source)
            .withListener(listener)
            .withContext(context)
            .withCode(code)
            .apply()
    }
}
