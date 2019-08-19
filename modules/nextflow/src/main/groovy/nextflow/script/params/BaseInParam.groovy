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

package nextflow.script.params

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.NF
import nextflow.exception.ProcessException
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.script.ProcessConfig
import nextflow.script.TokenVar
/**
 * Model a process generic input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
abstract class BaseInParam extends BaseParam implements InParam {

    protected fromObject

    protected bindObject

    protected owner

    /**
     * The channel to which the input value is bound
     */
    private inChannel

    /**
     * @return The input channel instance used by this parameter to receive the process inputs
     */
    DataflowReadChannel getInChannel() {
        init()
        return inChannel
    }

    BaseInParam(ProcessConfig config ) {
        this(config.getOwnerScript().getBinding(), config.getInputs())
    }

    /**
     * @param script The global script object
     * @param obj
     */
    BaseInParam( Binding binding, List holder, short ownerIndex = -1 ) {
        super(binding,holder,ownerIndex)
    }

    abstract String getTypeName()

    protected DataflowReadChannel inputValToChannel( value ) {
        checkFromNotNull(value)

        if( this instanceof DefaultInParam ) {
            assert value instanceof DataflowQueue
            return value
        }

        if ( value instanceof DataflowReadChannel || value instanceof DataflowBroadcast )  {
            return CH.getReadChannel(value)
        }

        if( NF.isDsl2() ) {
            final result = CH.value()
            result.bind(value)
            return result
        }

        // wrap any collections with a DataflowQueue
        if( value instanceof Collection ) {
            return CH.emitAndClose(CH.queue(), value)
        }

        // wrap any array with a DataflowQueue
        if ( value && value.class.isArray() ) {
            return CH.emitAndClose(CH.queue(), value as Collection)
        }

        value!=null ? CH.value(value) : CH.value() 
    }


    /**
     * Lazy parameter initializer.
     *
     * @return The parameter object itself
     */
    @Override
    protected void lazyInit() {

        if( fromObject == null && (bindObject == null || bindObject instanceof GString || bindObject instanceof Closure ) ) {
            throw new IllegalStateException("Missing 'bind' declaration in input parameter")
        }

        // fallback on the bind object if the 'fromObject' is not defined
        if( fromObject == null ) {
            fromObject = bindObject
        }

        // initialize the *inChannel* object based on the 'target' attribute
        def result
        if( fromObject instanceof TokenVar ) {
            // when the value is a variable reference
            // - use that name for the parameter itself
            // - get the variable value in the script binding
            result = getScriptVar(fromObject.name)
        }
        else if( fromObject instanceof Closure ) {
            result = fromObject.call()
        }
        else {
            result = fromObject
        }

        inChannel = inputValToChannel(result)
    }

    /**
     * @return The parameter name
     */
    String getName() {
        if( bindObject instanceof TokenVar )
            return bindObject.name

        if( bindObject instanceof String )
            return bindObject

        if( bindObject instanceof Closure )
            return '__$' + this.toString()

        throw new IllegalArgumentException("Invalid process input definition")
    }

    BaseInParam bind( def obj ) {
        this.bindObject = obj
        return this
    }

    private void checkFromNotNull(obj) {
        if( obj != null ) return
        def message = 'A process input channel evaluates to null'
        def name = null
        if( bindObject instanceof TokenVar )
            name = bindObject.name
        else if( bindObject instanceof CharSequence )
            name = bindObject.toString()
        if( name )
            message += " -- Invalid declaration `${getTypeName()} $name`"
        throw new IllegalArgumentException(message)
    }

    BaseInParam from( def obj ) {
        if(NF.isDsl2())
            throw new ScriptRuntimeException("Process clause `from` should not be provided when using DSL 2")
        setFrom(obj)
        return this
    }

    void setFrom( obj ) {
        checkFromNotNull(obj)
        fromObject = obj
    }

    Object getRawChannel() {
        if( CH.isChannel(fromObject) )
            return fromObject
        if( CH.isChannel(inChannel) )
            return inChannel
        throw new IllegalStateException("Missing input channel")
    }

    BaseInParam from( Object... obj ) {

        def normalize = obj.collect {
            if( it instanceof DataflowReadChannel )
                throw new IllegalArgumentException("Multiple channels are not allowed on 'from' input declaration")

            if( it instanceof Closure )
                return it.call()
            else
                it
        }

        fromObject = normalize as List
        return this
    }

    def decodeInputs( List inputs ) {
        final UNDEF = -1 as short
        def value = inputs[index]

        if( mapIndex == UNDEF || owner instanceof EachInParam )
            return value

        if( mapIndex != UNDEF ) {
            def result
            if( value instanceof Map ) {
                result = value.values()
            }
            else if( value instanceof Collection ) {
                result = value
            }
            else {
                result = [value]
            }

            try {
                return result[mapIndex]
            }
            catch( IndexOutOfBoundsException e ) {
                throw new ProcessException(e)
            }
        }

        return value
    }

}
