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
 */

package nextflow.script.params

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.extension.CH
import nextflow.script.ProcessConfig
import nextflow.script.TokenVar
import nextflow.util.ConfigHelper
/**
 * Model a process generic output parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseOutParam extends BaseParam implements OutParam {

    /** The out parameter name */
    protected String nameObj

    protected intoObj

    protected List<DataflowWriteChannel> outChannels = new ArrayList<>(10)

    @PackageScope
    boolean singleton

    String channelEmitName

    String channelTopicName

    BaseOutParam( Binding binding, List list, short ownerIndex = -1) {
        super(binding,list,ownerIndex)
    }

    BaseOutParam( ProcessConfig config ) {
        super(config.getOwnerScript().getBinding(), config.getOutputs())
    }

    Object clone() {
        final copy = (BaseOutParam)super.clone()
        copy.outChannels = new ArrayList<>(10)
        return copy
    }

    void lazyInit() {

        if( intoObj instanceof TokenVar || intoObj instanceof TokenVar[] ) {
            throw new IllegalArgumentException("Not a valid output channel argument: $intoObj")
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
        final channel = (target != null)
            ? outputValToChannel(target)
            : null

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
                channel = CH.create( singleton )

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


    BaseOutParam bind( def obj ) {
        if( obj instanceof TokenVar )
            this.nameObj = obj.name

        else
            this.nameObj = ( obj?.toString() ?: null )

        return this
    }

    void setInto( Object obj ) {
        intoObj = obj
    }

    DataflowWriteChannel getOutChannel() {
        init()
        return outChannels ? outChannels.get(0) : null
    }

    String getName() {
        if( nameObj != null )
            return nameObj.toString()
        throw new IllegalStateException("Missing 'name' property in output parameter")
    }

    @Override
    BaseOutParam setOptions(Map<String,?> opts) {
        super.setOptions(opts)
        return this
    }

    BaseOutParam setEmit( value ) {
        if( isNestedParam() )
            throw new IllegalArgumentException("Output `emit` option is not allowed in tuple components")
        if( !value )
            throw new IllegalArgumentException("Missing output `emit` name")
        if( !ConfigHelper.isValidIdentifier(value) ) {
            final msg = "Output emit '$value' is not a valid name -- Make sure it starts with an alphabetic or underscore character and it does not contain any blank, dot or other special characters"
            if( NF.strictMode )
                throw new IllegalArgumentException(msg)
            log.warn(msg)
        }
        this.channelEmitName = value
        return this
    }

    BaseOutParam setTopic( String name ) {
        if( isNestedParam() )
            throw new IllegalArgumentException("Output `topic` option it not allowed in tuple components")
        if( !name )
            throw new IllegalArgumentException("Missing output `topic` name")
        if( !ConfigHelper.isValidIdentifier(name) ) {
            final msg = "Output topic '$name' is not a valid name -- Make sure it starts with an alphabetic or underscore character and it does not contain any blank, dot or other special characters"
            throw new IllegalArgumentException(msg)
        }

        this.channelTopicName = name
        return this
    }
}
