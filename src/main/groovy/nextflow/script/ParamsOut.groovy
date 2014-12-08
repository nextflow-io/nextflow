/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.TaskConfig
/**
 * Model a process generic input parameter
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

interface OutParam {

    interface Mode {  }

    /**
     * @return The parameter name getter
     */
    String getName()

    /**
     * Defines the channel to which bind the output(s) in the script context
     *
     * @param value It can be a string representing a channel variable name in the script context. If
     *      the variable does not exist it creates a {@code DataflowVariable} in the script with that name.
     *      If the specified {@code value} is a {@code DataflowWriteChannel} object, use this object
     *      as the output channel
     * @return
     */
    OutParam into( def value )

    /**
     * @return The output channel instance
     */
    DataflowWriteChannel getOutChannel()

    short getIndex()

    Mode getMode()
}

enum BasicMode implements OutParam.Mode {

    standard, flatten;

    static OutParam.Mode parseValue( def x ) {

        if( x instanceof OutParam.Mode ) {
            return x
        }

        def value = x instanceof TokenVar ? x.name : ( x instanceof String ? x : null )
        if( value ) {
            return BasicMode.valueOf(value)
        }

        throw new IllegalArgumentException("Not a valid output 'mode' value: $value")

    }
}

@Slf4j
abstract class BaseOutParam extends BaseParam implements OutParam {

    /** The out parameter name */
    protected nameObj

    protected intoObj

    private outChannel

    protected OutParam.Mode mode = BasicMode.standard

    /** Whenever the channel has to closed on task termination */
    protected Boolean autoClose = Boolean.TRUE

    BaseOutParam( Binding binding, List list, short ownerIndex = -1) {
        super(binding,list,ownerIndex)
    }

    BaseOutParam( TaskConfig config ) {
        super(config.getOwnerScript().getBinding(), config.getOutputs())
    }

    void lazyInit() {
        def target
        if( intoObj instanceof TokenVar )
            target = intoObj.name

        else if( intoObj != null )
            target = intoObj

        else if( nameObj instanceof String )
            target = nameObj

        else
            throw new IllegalArgumentException("Missing 'into' in output param declaration")

        // define the output channel
        outChannel = outputValToChannel(target, DataflowQueue)
    }

    @groovy.transform.PackageScope
    BaseOutParam bind( def obj ) {
        if( obj instanceof TokenVar )
            this.nameObj = obj.name

        else
            this.nameObj = ( obj?.toString() ?: null )

        return this
    }

    BaseOutParam into( def value ) {
        intoObj = value
        return this
    }

    DataflowWriteChannel getOutChannel() {
        init()
        return outChannel
    }

    def String getName() {

        if( nameObj != null )
            return nameObj.toString()

        throw new IllegalStateException("Missing 'name' property in output parameter")
    }


    def BaseOutParam mode( def mode ) {
        this.mode = BasicMode.parseValue(mode)
        return this
    }

    OutParam.Mode getMode() { mode }

}


/**
 * Model a process *file* output parameter
 */
@Slf4j
@InheritConstructors
class FileOutParam extends BaseOutParam implements OutParam {

    /**
     * The character used to separate multiple names (pattern) in the output specification
     */
    protected String separatorChar = ':'

    /**
     * When {@code true} star wildcard (*) matches hidden files (files starting with a dot char)
     * By default it does not, coherently with linux bash rule
     */
    protected boolean includeHidden

    /**
     * When {@code true} file pattern includes input files as well as output files.
     * By default a file pattern matches only against files produced by the process, not
     * the ones received as input
     */
    protected boolean includeInputs

    private TokenGString gstring

    private Closure<String> dynamicObj

    String getSeparatorChar() { separatorChar }

    boolean getIncludeHidden() { includeHidden }

    boolean getIncludeInputs() { includeInputs }


    /**
     * @return {@code true} when the file name is parametric i.e contains a variable name to be resolved, {@code false} otherwise
     */
    boolean isParametric() { dynamicObj || gstring != null }

    FileOutParam separatorChar( String value ) {
        this.separatorChar = value
        return this
    }

    FileOutParam includeInputs( boolean flag ) {
        this.includeInputs = flag
        return this
    }

    FileOutParam includeHidden( boolean flag ) {
        this.includeHidden = flag
        return this
    }

    BaseOutParam bind( obj ) {

        if( obj instanceof TokenGString ) {
            log.warn "Parametric output file names should be defined with closures -- Replace `file \"${obj.text}\"` with `file { \"${obj.text}\" }`"
            gstring = obj
            return this
        }

        if( obj instanceof TokenVar ) {
            this.nameObj = obj.name
            dynamicObj = { delegate.containsKey(obj.name) ? delegate.get(obj.name): obj.name }
            return this
        }

        if( obj instanceof Closure ) {
            dynamicObj = obj
            return this
        }

        // fallback on super class
        super.bind(obj)
    }

    List<String> getFilePatterns(Map context) {

        def nameString
        if( dynamicObj ) {
            nameString = context.with(dynamicObj)
        }
        else if( gstring ) {
            def strict = getName() == null
            nameString = gstring.resolve { String it -> resolveName(context, it, strict) }
        }
        else {
            nameString = nameObj
        }

        return separatorChar ? nameString.split(/\${separatorChar}/) : [nameString]
    }

    /**
     * Override the default to allow null as a value name
     * @return
     */
    String getName() {
        return nameObj ? super.getName() : null
    }

}


/**
 * Model a process *value* output parameter
 */
@InheritConstructors
class ValueOutParam extends BaseOutParam { }

/**
 * Model the process *stdout* parameter
 */
@InheritConstructors
class StdOutParam extends BaseOutParam { }


@InheritConstructors
class SetOutParam extends BaseOutParam {

    enum CombineMode implements OutParam.Mode { combine }

    final List<BaseOutParam> inner = []

    String getName() { toString() }

    SetOutParam bind( Object... obj ) {

        obj.each { item ->
            if( item instanceof TokenVar )
                create(ValueOutParam).bind(item)

            else if( item instanceof TokenValCall )
                create(ValueOutParam).bind(item.name)

            else if( item instanceof TokenStdoutCall || item == '-'  )
                create(StdOutParam).bind('-')

            else if( item instanceof String )
                create(FileOutParam).bind(item)

            else if( item instanceof TokenFileCall )
                // note that 'filePattern' can be a string or a Gstring
                create(FileOutParam).bind(item.target)

            else if( item instanceof TokenGString )
                create(FileOutParam).bind(item)

            else
                throw new IllegalArgumentException("Invalid map output parameter declaration -- item: ${item}")
        }

        return this
    }

    def void lazyInit() {
        if( intoObj == null )
            throw new IllegalStateException("Missing 'into' channel in output parameter declaration")

        super.lazyInit()
    }

    protected <T extends BaseOutParam> T create(Class<T> type) {
        type.newInstance(binding,inner,index)
    }

    def SetOutParam mode( def value ) {

        def str = value instanceof String ? value : ( value instanceof TokenVar ? value.name : null )
        if( str ) {
            try {
                this.mode = CombineMode.valueOf(str)
            }
            catch( Exception e ) {
                super.mode(value)
            }
        }

        return this
    }
}

/**
 * Container to hold all process outputs
 */
class OutputsList implements List<OutParam> {

    @Delegate
    private List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() {
        target.collect { OutParam it -> it.getOutChannel() }
    }

    List<String> getNames() { target *. name }

    def <T extends OutParam> List<T> ofType( Class<T>... classes ) {
        (List<T>) target.findAll { it.class in classes }
    }

}
