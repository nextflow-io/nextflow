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
import java.nio.file.Path

import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.exception.IllegalFileException
import nextflow.file.FilePatternSplitter
import nextflow.processor.ProcessConfig
import nextflow.util.BlankSeparatedList
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

    List<DataflowWriteChannel> getOutChannels()

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

    private List<DataflowWriteChannel> outChannels = []

    protected OutParam.Mode mode = BasicMode.standard

    @PackageScope
    boolean singleton

    BaseOutParam( Binding binding, List list, short ownerIndex = -1) {
        super(binding,list,ownerIndex)
    }

    BaseOutParam( ProcessConfig config ) {
        super(config.getOwnerScript().getBinding(), config.getOutputs())
    }

    void lazyInit() {

        if( intoObj instanceof TokenVar[] ) {
            intoObj.each { lazyInitImpl(it) }
        }
        else if( intoObj != null ) {
            lazyInitImpl(intoObj)
        }
        else if( nameObj instanceof String ) {
            lazyInitImpl(nameObj)
        }

    }

    @PackageScope
    void setSingleton( boolean value ) {
        this.singleton = value
    }

    @PackageScope
    void lazyInitImpl( def target ) {
        def channel = null
        if( target instanceof TokenVar ) {
            channel = outputValToChannel(target.name)
        }
        else if( target != null ) {
            channel = outputValToChannel(target)
        }

        if( channel ) {
            outChannels.add(channel)
        }
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
    final protected DataflowWriteChannel outputValToChannel( Object channel ) {

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
                    log.trace "Creating new output channel > $local"
                }
                else {
                    log.warn "Output channel `$local` overrides another variable with the same name declared in the script context -- Rename it to avoid possible conflicts"
                }

                // instantiate the new channel
                channel = singleton && mode==BasicMode.standard ? new DataflowVariable() : new DataflowQueue()

                // bind it to the script on-fly
                if( local != '-' && binding ) {
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

    BaseOutParam into( TokenVar... vars ) {
        intoObj = vars
        return this
    }

    @Deprecated
    DataflowWriteChannel getOutChannel() {
        init()
        return outChannels ? outChannels.get(0) : null
    }

    List<DataflowWriteChannel> getOutChannels() {
        init()
        return outChannels
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
 * Implements an optional file output option
 */
trait OptionalParam {

    private boolean optional

    boolean getOptional() { optional }

    def optional( boolean value ) {
        this.optional = value
        return this
    }

}

/**
 * Placeholder trait to mark a missing optional output parameter
 */
trait MissingParam {

    OutParam missing

}

/**
 * Model a process *file* output parameter
 */
@Slf4j
@InheritConstructors
class FileOutParam extends BaseOutParam implements OutParam, OptionalParam {

    /**
     * ONLY FOR TESTING DO NOT USE
     */
    protected FileOutParam(Map params) {
        super(new Binding(), [])
    }

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

    protected boolean glob = true

    private GString gstring

    private Closure<String> dynamicObj

    private String filePattern

    String getSeparatorChar() { separatorChar }

    boolean getHidden() { includeHidden }

    boolean getIncludeInputs() { includeInputs }

    String getType() { type }

    Integer getMaxDepth() { maxDepth }

    boolean getFollowLinks() { followLinks }

    boolean getGlob() { glob }


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

    FileOutParam glob( boolean value ) {
        glob = value
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

        this.filePattern = obj.toString()
        return this
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
            entry = filePattern
        }

        if( !entry )
            return []

        if( entry instanceof Path )
            return [ relativize(entry, workDir) ]

        // handle a collection of files
        if( entry instanceof BlankSeparatedList || entry instanceof List ) {
            return entry.collect { relativize(it.toString(), workDir) }
        }

        // normalize to a string object
        final nameString = entry.toString()
        if( separatorChar && nameString.contains(separatorChar) ) {
            return nameString.split(/\${separatorChar}/).collect { String it-> relativize(it, workDir) }
        }

        return [relativize(nameString, workDir)]

    }

    @PackageScope String getFilePattern() { filePattern }

    @PackageScope
    static String clean(String path) {
        while (path.startsWith('/') ) {
            path = path.substring(1)
        }
        return path
    }

    @PackageScope
    String relativize(String path, Path workDir) {
        if( !path.startsWith('/') )
            return path

        final dir = workDir.toString()
        if( !path.startsWith(dir) )
            throw new IllegalFileException("File `$path` is out of the scope of process working dir: $workDir")

        if( path.length()-dir.length()<2 )
            throw new IllegalFileException("Missing output file name")

        return path.substring(dir.size()+1)
    }

    @PackageScope
    String relativize(Path path, Path workDir) {
        if( !path.isAbsolute() )
            return glob ? FilePatternSplitter.GLOB.escape(path) : path

        if( !path.startsWith(workDir) )
            throw new IllegalFileException("File `$path` is out of the scope of process working dir: $workDir")

        if( path.nameCount == workDir.nameCount )
            throw new IllegalFileException("Missing output file name")

        final rel = path.subpath(workDir.getNameCount(), path.getNameCount())
        return glob ? FilePatternSplitter.GLOB.escape(rel) : rel
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
class SetOutParam extends BaseOutParam implements OptionalParam {

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
                // note that 'filePattern' can be a string or a GString
                create(FileOutParam).bind(item.target)

            else
                throw new IllegalArgumentException("Invalid `set` output parameter declaration -- item: ${item}")
        }

        return this
    }

    protected <T extends BaseOutParam> T create(Class<T> type) {
        type.newInstance(binding,inner,index)
    }

    @Override
    void lazyInit() {
        super.lazyInit()
        inner.each { opt ->
            if( opt instanceof FileOutParam ) opt.optional(this.optional)
        }
    }

    SetOutParam mode( def value ) {

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

final class DefaultOutParam extends StdOutParam {

    DefaultOutParam( ProcessConfig config ) {
        super(config)
        bind('-')
        into(new DataflowQueue())
    }
}

/**
 * Container to hold all process outputs
 */
class OutputsList implements List<OutParam> {

    @Delegate
    private List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() {
        final List<DataflowWriteChannel> result = []
        target.each { OutParam it -> result.addAll(it.getOutChannels()) }
        return result
    }

    List<String> getNames() { target *. name }

    def <T extends OutParam> List<T> ofType( Class<T>... classes ) {
        (List<T>) target.findAll { it.class in classes }
    }

    void setSingleton( boolean value ) {
        target.each { BaseOutParam param -> param.singleton = value }
    }
}
