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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.script.ProcessConfig
import nextflow.script.TokenVar

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

    protected OutParam.Mode mode = BasicMode.standard

    @PackageScope
    boolean singleton

    String channelEmitName

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

        if( intoObj instanceof TokenVar[] ) {
            if( NF.dsl2 )
                throw new IllegalArgumentException("Not a valid output channel argument: $intoObj")
            for( def it : intoObj ) { lazyInitImpl(it) }
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
            assert !NF.dsl2
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
                channel = CH.create( singleton && mode==BasicMode.standard )

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

    BaseOutParam into( def value ) {
        if( NF.dsl2 )
            throw new ScriptRuntimeException("Process clause `into` should not be provided when using DSL 2")
        this.intoObj = value
        return this
    }

    BaseOutParam into( TokenVar... vars ) {
        if( NF.dsl2 )
            throw new ScriptRuntimeException("Process clause `into` should not be provided when using DSL 2")
        intoObj = vars
        return this
    }

    void setInto( Object obj ) {
        intoObj = obj
    }

    DataflowWriteChannel getOutChannel() {
        init()
        return outChannels ? outChannels.get(0) : null
    }

    @Deprecated
    List<DataflowWriteChannel> getOutChannels() {
        init()
        return outChannels
    }

    String getName() {
        if( nameObj != null )
            return nameObj.toString()
        throw new IllegalStateException("Missing 'name' property in output parameter")
    }


    BaseOutParam mode( def mode ) {
        if( NF.isDsl2() )
            log.warn "Process output `mode` is not supported any more"
        this.mode = BasicMode.parseValue(mode)
        return this
    }

    OutParam.Mode getMode() { mode }

    @Override
    BaseOutParam setOptions(Map<String,?> opts) {
        super.setOptions(opts)
        return this
    }

    BaseOutParam setEmit( value ) {
        if( isNestedParam() )
            throw new IllegalArgumentException("Output `emit` option it not allowed in tuple components")
        this.channelEmitName = value
        return this
    }
}
