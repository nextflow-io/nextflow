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

package nextflow.script

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Nextflow
import nextflow.exception.ProcessException
import nextflow.extension.ToListOp
import nextflow.processor.ProcessConfig
/**
 * Base class for input/output parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseParam {

    final protected Binding binding

    final protected List<BaseParam> holder

    final short index

    final short mapIndex

    private boolean initialized

    BaseParam ( Binding binding, List holder, int ownerIndex = -1 ) {
        this.binding = binding
        this.holder = holder

        /*
         * by default the index is got from 'holder' current size
         * and the mapIndex is =1 (not defined)
         */
        if( ownerIndex == -1 ) {
            index = holder.size()
            mapIndex = -1
        }

        /*
         * when the owner index is provided (not -1) it is used as
         * the main index and the map index is got from the 'holder' size
         */
        else {
            index = ownerIndex
            mapIndex = holder.size()
        }

        // add the the param to the holder list
        holder.add(this)
    }

    String toString() {
        def p = mapIndex == -1 ? index : "$index:$mapIndex"
        return "${this.class.simpleName.toLowerCase()}<$p>"
    }

    /**
     * Lazy initializer
     */
    protected abstract void lazyInit()

    /**
     * Initialize the parameter fields if needed
     */
    final protected void init() {
        if( initialized ) return
        lazyInit()

        // flag as initialized
        initialized = true
    }


    /**
     * Get the value of variable {@code name} in the script context
     *
     * @param name The variable name
     * @param strict If {@code true} raises a {@code MissingPropertyException} when the specified variable does not exist
     * @return The variable object
     */
    protected getScriptVar(String name, boolean strict ) {
        if( binding.hasVariable(name) ) {
            return binding.getVariable(name)
        }

        if( strict )
            throw new MissingPropertyException(name,this.class)

        return null
    }

    protected getScriptVar( String name ) {
        getScriptVar(name,true)
    }

}

/**
 * Basic interface for *all* input parameters
 */
interface InParam {

    String getName()

    DataflowReadChannel getInChannel()

    InParam from( Object value )

    InParam from( Object... values )

    short index

    short mapIndex

    def decodeInputs( List values )

}

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

    BaseInParam( ProcessConfig config ) {
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

        if ( value instanceof DataflowBroadcast )  {
            return value.createReadChannel()
        }

        if( value instanceof DataflowReadChannel ) {
            return value
        }

        // wrap any collections with a DataflowQueue
        if( value instanceof Collection ) {
            return Nextflow.channel(value as List)
        }

        // wrap any array with a DataflowQueue
        if ( value && value.class.isArray() ) {
            return Nextflow.channel(value as List)
        }

        // wrap a single value with a DataflowVariable
        return Nextflow.variable(value)

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
        checkFromNotNull(obj)
        fromObject = obj
        return this
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

/**
 *  Represents a process *file* input parameter
 */
@Slf4j
@InheritConstructors
class FileInParam extends BaseInParam  {

    protected filePattern

    @Override String getTypeName() { 'file' }

    /**
     * Define the file name
     */
    FileInParam name( obj ) {
        if( obj instanceof String ) {
            filePattern = obj
            return this
        }

        if( obj instanceof GString ) {
            filePattern = obj
            return this
        }

        // the ability to pass a closure as file name has been replaced by
        // lazy gstring -- this should be deprecated
        if( obj instanceof Closure ) {
            filePattern = obj
            return this
        }

        throw new IllegalArgumentException()
    }

    String getName() {

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return entry?.key
        }

        if( bindObject instanceof GString ) {
            return '__$' + this.toString()
        }

        return super.getName()

    }

    String getFilePattern(Map ctx = null) {

        if( filePattern != null  )
            return resolve(ctx,filePattern)

        if( bindObject instanceof Map ) {
            def entry = bindObject.entrySet().first()
            return resolve(ctx, entry?.value)
        }

        if( bindObject instanceof TokenVar )
            return filePattern = '*'

        if( bindObject != null )
            return resolve(ctx, bindObject)

        return filePattern = '*'
    }

    private resolve( Map ctx, value ) {
        if( value instanceof GString ) {
            value.cloneWith(ctx)
        }

        else if( value instanceof Closure ) {
            return ctx.with(value)
        }

        else
            return value
    }


}

/**
 *  Represents a process *environment* input parameter
 */
@InheritConstructors
class EnvInParam extends BaseInParam {
    @Override String getTypeName() { 'env' }
}

/**
 *  Represents a process *value* input parameter
 */
@InheritConstructors
class ValueInParam extends BaseInParam {
    @Override String getTypeName() { 'val' }
}

/**
 *  Represents a process *stdin* input parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class StdInParam extends BaseInParam {

    String getName() { '-' }

    @Override String getTypeName() { 'stdin' }
}

/**
 *  Represents a process input *iterator* parameter
 */
@InheritConstructors
@Slf4j
class EachInParam extends BaseInParam {

    @Override String getTypeName() { 'each' }

    private List<InParam> inner = []

    String getName() { '__$'+this.toString() }

    EachInParam bind( def obj ) {
        def nested = ( obj instanceof TokenFileCall
                    ? new FileInParam(binding, inner, index).bind(obj.target)
                    : new ValueInParam(binding, inner, index).bind(obj) )
        nested.owner = this
        this.bindObject = nested.bindObject
        return this
    }

    InParam getInner() { inner[0] }

    @Override
    protected DataflowReadChannel inputValToChannel( value ) {
        def variable = normalizeToVariable( value )
        super.inputValToChannel(variable)
    }

    @PackageScope
    DataflowReadChannel normalizeToVariable( value ) {
        def result
        if( value instanceof DataflowExpression ) {
            result = value
        }

        else if( value instanceof DataflowReadChannel ) {
            result = ToListOp.apply(value)
        }

        else {
            result = new DataflowVariable()
            result.bind(value)
        }

        return result.chainWith { it instanceof Collection || it == null ? it : [it] }
    }

}

@InheritConstructors
class SetInParam extends BaseInParam {

    final List<InParam> inner = []

    @Override String getTypeName() { 'set' }

    String getName() { '__$'+this.toString() }

    SetInParam bind( Object... obj ) {

        obj.each { item ->

            if( item instanceof TokenVar )
                newItem(ValueInParam).bind(item)

            else if( item instanceof TokenFileCall )
                newItem(FileInParam).bind( item.target )

            else if( item instanceof Map )
                newItem(FileInParam).bind(item)

            else if( item instanceof TokenValCall )
                newItem(ValueInParam).bind(item.val)

            else if( item instanceof TokenEnvCall )
                newItem(EnvInParam).bind(item.val)

            else if( item instanceof TokenStdinCall )
                newItem(StdInParam)

            else if( item instanceof GString )
                newItem(FileInParam).bind(item)

            else if( item == '-' )
                newItem(StdInParam)

            else if( item instanceof String )
                newItem(FileInParam).bind(item)

            else
                throw new IllegalArgumentException()
        }

        return this

    }

    private <T extends BaseInParam> T newItem( Class<T> type )  {
        type.newInstance(binding, inner, index)
    }

}

final class DefaultInParam extends ValueInParam {

    @Override String getTypeName() { 'default' }

    DefaultInParam(ProcessConfig config) {
        super(config)
        final channel = new DataflowQueue(); channel.bind(Boolean.TRUE)
        from(channel)
        bind('$')
    }
}


/**
 * Container to hold all process outputs
 */
@Slf4j
class InputsList implements List<InParam> {

    @Delegate
    private List<InParam> target = new LinkedList<>()

    List<DataflowReadChannel> getChannels() {
        target.collect { InParam it -> it.getInChannel() }
    }

    List<String> getNames() { target *. name }


    def <T extends InParam> List<T> ofType( Class<T> clazz ) {
        (List<T>) target.findAll { it.class == clazz }
    }

    boolean allScalarInputs() {
        for( InParam param : target ) {
            if( param.inChannel instanceof DataflowQueue )
                return false
        }
        return true
    }

}

