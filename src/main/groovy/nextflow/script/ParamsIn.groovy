/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.script

import java.util.concurrent.atomic.AtomicBoolean

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Nextflow
import nextflow.extension.DataflowExtensions
import nextflow.processor.ProcessConfig
/**
 * Base class for input/output parameters
 *
 * @author Paolo DI Tommaso <paolo.ditommaso@gmail.com>
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


    protected DataflowReadChannel inputValToChannel( def value ) {

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
     * Creates a channel variable in the script context
     *
     * @param channel it can be a string representing a channel variable name in the script context. If
     *      the variable does not exist it creates a {@code DataflowVariable} in the script with that name.
     *      If the specified {@code value} is a {@code DataflowWriteChannel} object, use this object
     *      as the output channel
     *
     * @param factory The type of the channel to create, either {@code DataflowVariable} or {@code DataflowQueue}
     * @return The created (or specified) channel instance
     */
    final protected DataflowWriteChannel outputValToChannel( Object channel, Class<DataflowWriteChannel> factory ) {

        if( channel instanceof String ) {
            // the channel is specified by name
            def local = channel

            // look for that name in the 'script' context
            channel = binding.hasVariable(local) ? binding.getVariable(local) : null
            if( channel instanceof DataflowWriteChannel ) {
                // that's OK -- nothing to do
            }
            else {
                if( channel == null ) {
                    log.debug "output > channel unknown: $local -- creating a new instance"
                }
                else {
                    log.warn "Duplicate output channel name: '$channel' in the script context -- it's worth to rename it to avoid possible conflicts"
                }

                // instantiate the new channel
                channel = factory.newInstance()

                // bind it to the script on-fly
                if( local != '-' && binding) {
                    // bind the outputs to the script scope
                    binding.setVariable(local, channel)
                }
            }
        }

        if( channel instanceof DataflowWriteChannel ) {
            return channel
        }

        throw new IllegalArgumentException("Invalid output channel reference")
    }

    @Deprecated
    final protected resolveName( Map context, String name, boolean strict = true ) {
        if( context && context.containsKey(name) )
            return context.get(name)

//        if( binding.hasVariable(name) )
//            return binding.getVariable(name)

        if( strict )
            throw new MissingPropertyException(name)

        return name
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
    def String getName() {
        if( bindObject instanceof TokenVar )
            return bindObject.name

        if( bindObject instanceof String )
            return bindObject

        if( bindObject instanceof Closure )
            return '__$' + this.toString()

        throw new IllegalArgumentException()
    }

    BaseInParam bind( def obj ) {
        this.bindObject = obj
        return this
    }

    BaseInParam from( def obj ) {
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



}

/**
 *  Represents a process *file* input parameter
 */
@Slf4j
@InheritConstructors
class FileInParam extends BaseInParam  {

    protected filePattern

    private AtomicBoolean warnShown = new AtomicBoolean()

    private void warn(String str) {
        if(!warnShown.getAndSet(true))
        log.warn "Dynamic input file name has to be defined using a closure -- Replace `file \"${str}\"` with `file { \"${str}\" }`"
    }

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


}

/**
 *  Represents a process *value* input parameter
 */
@InheritConstructors
class ValueInParam extends BaseInParam { }

/**
 *  Represents a process *stdin* input parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class StdInParam extends BaseInParam {

    String getName() { '-' }

}

/**
 *  Represents a process input *iterator* parameter
 */
@InheritConstructors
@Slf4j
class EachInParam extends BaseInParam {

    @Override
    protected DataflowReadChannel inputValToChannel( value ) {

        def variable = normalizeToVariable( value, name )
        super.inputValToChannel(variable)
    }

    @PackageScope
    static DataflowVariable normalizeToVariable( value, String name = null) {
        if( value instanceof DataflowReadChannel ) {
            log.warn "Using queue channel on each parameter declaration should be avoided -- take in consideration to change declaration for each: '$name' parameter"
            // everything is mapped to a collection
            value = readValue(value)
        }

        if( !(value instanceof Collection) ) {
            value = [value]
        }

        // the collection is wrapped to a "scalar" dataflow variable
        Nextflow.variable(value)
    }


    /**
     * Converts the specified arguments to a {@code List} data type
     *
     * @param item
     * @return
     */
    @PackageScope
    static readValue( DataflowReadChannel channel ) {

        if( channel instanceof DataflowExpression ) {
            return channel.getVal()
        }
        else {
            return DataflowExtensions.toList(channel).getVal()
        }

    }

}

@InheritConstructors
class SetInParam extends BaseInParam {

    final List<InParam> inner = []

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
                newItem(ValueInParam).bind(item.name)

            else if( item instanceof TokenEnvCall )
                newItem(EnvInParam).bind(item.name)

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

