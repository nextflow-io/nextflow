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

import java.nio.file.Path

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.IllegalFileException
import nextflow.processor.ProcessConfig
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

    BaseOutParam( ProcessConfig config ) {
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

    @PackageScope
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
     * ONLY FOR TESTING DO NOT USE
     */
    protected FileOutParam(Map params) { }

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

    /**
     * The type of path to output, either {@code file}, {@code dir} or {@code any}
     */
    protected String type

    /**
     * Maximum number of directory levels to visit (default: no limit)
     */
    protected Integer maxDepth

    /**
     * When true it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: true)
     */
    protected boolean followLinks = true

    private GString gstring

    private Closure<String> dynamicObj

    String getSeparatorChar() { separatorChar }

    boolean getHidden() { includeHidden }

    boolean getIncludeInputs() { includeInputs }

    String getType() { type }

    Integer getMaxDepth() { maxDepth }

    boolean getFollowLinks() { followLinks }


    /**
     * @return {@code true} when the file name is parametric i.e contains a variable name to be resolved, {@code false} otherwise
     */
    boolean isDynamic() { dynamicObj || gstring != null }

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

    FileOutParam hidden( boolean flag ) {
        this.includeHidden = flag
        return this
    }

    FileOutParam type( String value ) {
        assert value in ['file','dir','any']
        type = value
        return this
    }

    FileOutParam maxDepth( int value ) {
        maxDepth = value
        return this
    }

    FileOutParam followLinks( boolean value ) {
        followLinks = value
        return this
    }


    BaseOutParam bind( obj ) {

        if( obj instanceof GString ) {
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

    List<String> getFilePatterns(Map context, Path workDir) {

        def entry = null
        if( dynamicObj ) {
            entry = context.with(dynamicObj)
        }
        else if( gstring != null ) {
            def strict = (getName() == null)
            try {
                entry = gstring.cloneWith(context)
            }
            catch( MissingPropertyException e ) {
                if( strict )
                    throw e
            }
        }
        else {
            entry = nameObj
        }

        if( !entry )
            return []

        if( entry instanceof Path )
            return [ relativizeName(entry, workDir) ]

        // normalize to a string object
        final nameString = entry.toString()
        if( separatorChar && nameString.contains(separatorChar) )
            return nameString.split(/\${separatorChar}/).collect { String it-> clean(it) }

        return [clean(nameString)]

    }

    @PackageScope
    static String clean(String path) {
        while (path.startsWith('/') ) {
            path = path.substring(1)
        }
        return path
    }

    @PackageScope
    static String relativizeName(Path path, Path workDir) {
        if( path.isAbsolute() ) {
            if( !path.startsWith(workDir) ) {
                throw new IllegalFileException("File `$path` is out of the scope of process working dir: $workDir")
            }

            return path.subpath(workDir.getNameCount(), path.getNameCount()).toString()
        }

        return path.toString()
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
class ValueOutParam extends BaseOutParam {

    protected target

    String getName() {
        return nameObj ? super.getName() : null
    }


    BaseOutParam bind( def obj ) {
        // the target value object
        target = obj

        // retrieve the variable name to be used to fetch the value
        if( obj instanceof TokenVar ) {
            this.nameObj = obj.name
        }

        return this
    }

    /**
     * Given the {@link nextflow.processor.TaskContext} object resolve the actual value
     * to which this param is bound
     *
     * @param context An instance of {@link nextflow.processor.TaskContext} holding the task evaluation context
     * @return The actual value to which this out param is bound
     */
    def resolve( Map context ) {

        switch( target ) {
        case TokenVar:
            return context.get(target.name)

        case Closure:
            return target.cloneWith(context).call()

        case GString:
            return target.cloneWith(context).toString()

        default:
            return target
        }
    }

}

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
                create(ValueOutParam).bind(item.val)

            else if( item instanceof GString )
                create(FileOutParam).bind(item)

            else if( item instanceof TokenStdoutCall || item == '-'  )
                create(StdOutParam).bind('-')

            else if( item instanceof String )
                create(FileOutParam).bind(item)

            else if( item instanceof TokenFileCall )
                // note that 'filePattern' can be a string or a Gstring
                create(FileOutParam).bind(item.target)

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
