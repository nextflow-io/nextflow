/*
 * Copyright (c) 2012, the authors.
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

package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
/**
 * Model a process generic input parameter
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

interface OutParam {

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
    OutParam to( def value )

    /**
     * @return The output channel instance
     */
    DataflowWriteChannel getOutChannel()

}

@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class BaseOutParam extends BaseParam implements OutParam {

    /** The out parameter name */
    protected String name

    protected Object outTarget

    private outChannel

    /** Whenever the channel has to closed on task termination */
    protected Boolean autoClose = Boolean.TRUE

    BaseOutParam( Script script, String name ) {
        super(script)
        this.name = name
    }

    BaseOutParam lazyInit() {
        def value = outTarget ?: name
        outChannel = outputValToChannel(value, DataflowQueue)
        return this
    }

    BaseOutParam to( def value ) {
        this.outTarget = value
        return this
    }

    DataflowWriteChannel getOutChannel() {
        init()
        return outChannel
    }

    OutParam autoClose( boolean value ) {
        this.autoClose = value
        return this
    }

    def String getName() { name }

    def Boolean getAutoClose() { autoClose }

}


/**
 * Model a process *file* output parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true, includeNames = true)
class FileOutParam extends BaseOutParam implements OutParam {

    /**
     * Whenever multiple files matching the same name pattern have to be output as grouped collection
     * or bound as single entries
     */
    protected boolean flat = false

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

    boolean getFlat() { flat }

    String getSeparatorChar() { separatorChar }

    boolean getIncludeHidden() { includeHidden }

    boolean getIncludeInputs() { includeInputs }

    FileOutParam flat( boolean value ) {
        this.flat = value
        return this
    }

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
}


/**
 * Model a process *value* output parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true, includeNames = true)
class ValueOutParam extends BaseOutParam { }

/**
 * Model the process *stdout* parameter
 */
@ToString(includePackage=false, includeSuper = true, includeNames = true)
class StdOutParam extends BaseOutParam {
    StdOutParam(Script script) { super(script,'-') }
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

    def <T extends OutParam> List<T> ofType( Class<T> clazz ) { (List<T>) target.findAll { it.class == clazz } }

    void eachParam (Closure closure) {
        target.each { OutParam param -> closure.call(param.name, param.getOutChannel()) }
    }

}
