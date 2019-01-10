/*
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

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.operator.PoisonPill
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import org.codehaus.groovy.runtime.callsite.BooleanReturningMethodInvoker
import static nextflow.extension.DataflowHelper.newOperator
import static nextflow.util.CheckHelper.checkParams
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BufferOp {

    static final private Map BUFFER_PARAMS = [
            skip: Integer,
            size: Integer,
            remainder: Boolean
    ]

    @Lazy
    private static Session session = Global.session as Session

    private Map params

    private DataflowReadChannel source

    private DataflowQueue target

    private Object openCondition

    private Object closeCondition

    private int skip

    private int size

    private boolean remainder

    private int skipCount = 0

    private int itemCount = 0


    BufferOp( DataflowReadChannel source ) {
        assert source
        this.source = source
    }

    BufferOp setParams( Map params ) {
        checkParams( 'buffer', params, BUFFER_PARAMS )
        this.params = params
        return this
    }

    BufferOp setStartCriteria( Object criteria ) {
        this.openCondition = criteria
        return this
    }

    BufferOp setCloseCriteria( Object criteria ) {
        this.closeCondition = criteria
        return this
    }

    DataflowQueue apply() {
        target = new DataflowQueue()

        if( params?.skip )
            this.skip = params.skip as int
        if( params?.size )
            this.size = params.size as int 
        if( params?.remainder == true )
            this.remainder = true

        if( (skip||size) && openCondition )
            throw new IllegalArgumentException()

        if( (skip||size) && closeCondition )
            throw new IllegalArgumentException()

        Closure c1=null
        Closure c2=null
        if( skip || size ) {
            c1 = createSkipCriteria()
            c2 = createSizeCriteria()
        }
        else {
            if( openCondition ) c1 = createCriteria(openCondition)
            if( closeCondition ) c2 = createCriteria(closeCondition)
        }

        buffer0(source, target, c1, c2, remainder)
        return target
    }

    private Closure createCriteria( Object condition ) {
        def invoker = new BooleanReturningMethodInvoker("isCase")
        return { Object it -> invoker.invoke(condition, it) }
    }

    private Closure createSkipCriteria( ) {
        return {
            skipCount +=1
            if( skipCount > skip ) {
                skipCount = 0
                return true
            }
            return false
        }
    }

    private Closure createSizeCriteria( ) {
        return {
            itemCount +=1
            if( itemCount-skip == size ) {
                itemCount = 0;
                return true
            }
            return false
        }
    }

    @CompileDynamic
    static private <V> void buffer0(DataflowReadChannel<V> source, DataflowQueue target, Closure startingCriteria, Closure closeCriteria, boolean remainder ) {
        assert closeCriteria

        // the list holding temporary collected elements
        def buffer = []
        final stopOnFirst = source instanceof DataflowExpression

        // -- intercepts the PoisonPill and sent out the items remaining in the buffer when the 'remainder' flag is true
        final listener = new DataflowEventAdapter() {

            @Override
            Object controlMessageArrived(final DataflowProcessor processor, final DataflowReadChannel<Object> channel, final int index, final Object message) {
                if( message instanceof PoisonPill && remainder && buffer.size() ) {
                    target.bind(buffer)
                }
                return message
            }

            @Override
            void afterRun(DataflowProcessor processor, List<Object> messages) {
                if( !stopOnFirst )
                    return
                if( remainder && buffer)
                    target.bind(buffer)
                target.bind(Channel.STOP)
            }

            @Override
            boolean onException(DataflowProcessor processor, Throwable e) {
                DataflowExtensions.log.error("@unknown", e)
                session.abort(e)
                return true
            }
        }

        // -- open frame flag
        boolean isOpen = startingCriteria == null

        // -- the operator collecting the elements
        newOperator( source, target, listener ) {
            if( isOpen ) {
                buffer << it
            }
            else if( startingCriteria.call(it) ) {
                isOpen = true
                buffer << it
            }

            if( closeCriteria.call(it) ) {
                ((DataflowProcessor) getDelegate()).bindOutput(buffer);
                buffer = []
                // when a *startingCriteria* is defined, close the open frame flag
                isOpen = (startingCriteria == null)
            }

        }
    }

}
