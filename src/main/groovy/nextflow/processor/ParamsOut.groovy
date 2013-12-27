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

@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class OutParam {

    final Script script

    /** The out parameter name */
    protected String name

    protected Object using

    /** Whenever the channel has to closed on task termination */
    protected Boolean autoClose = Boolean.TRUE


    /** The channel over which entries are sent */
    @Lazy
    DataflowWriteChannel channel = {

        return using ? channelRef(using) : channelRef(name?.replaceAll(/\./,'_'))

    } ()



    OutParam( Script script, String name ) {
        this.name = name
        this.script = script
    }

    OutParam using( def value ) {
        this.using = value
        return this
    }

    protected DataflowWriteChannel channelRef( Object channel ) {
        log.trace "output using > channel ref: $channel"
        channelRef(script, channel, { new DataflowQueue() })
    }

    @groovy.transform.PackageScope
    static DataflowWriteChannel channelRef( Script script, Object channel, Closure<DataflowWriteChannel> factory ) {

        if( channel instanceof String ) {
            // the channel is specified by name
            def local = channel

            def binding = script.getBinding()

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
                channel = factory.call()

                // bind it to the script on-fly
                if( local != '-' && script) {
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
class FileOutParam extends OutParam {

    /**
     * Whenever multiple files matching the same name pattern have to be output as grouped collection
     * or bound as single entries
     */
    protected Boolean joint = Boolean.FALSE

    /**
     * The character used to separate multiple names (pattern) in the output specification
     */
    protected String separatorChar = ':'

    protected boolean includeHidden

    protected boolean includeInputs

    boolean getJoint() { joint }

    String getSeparatorChar() { separatorChar }

    boolean getIncludeHidden() { includeHidden }

    boolean getIncludeInputs() { includeInputs }

    FileOutParam joint( boolean value ) {
        this.joint = value
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
class ValueOutParam extends OutParam { }

/**
 * Model the process *stdout* parameter
 */
@ToString(includePackage=false, includeSuper = true, includeNames = true)
class StdOutParam extends OutParam { StdOutParam(Script script) { super(script,'-') } }

/**
 * Container to hold all process outputs
 */
class OutputsList implements List<OutParam> {

    @Delegate
    List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() { target *.channel }

    List<String> getNames() { target *. name }

    def <T extends OutParam> List<T> ofType( Class<T> clazz ) { (List<T>) target.findAll { it.class == clazz } }

    void eachParam (Closure closure) {
        target.each { OutParam param -> closure.call(param.name, param.channel) }
    }

}
