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
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Nextflow
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Holds the task configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskConfig extends GroovyObjectSupport {

    private Map configProperties = new LinkedHashMap()

    /**
     * Initialize the config object with the defaults values
     */
    TaskConfig() {

        configProperties.with {
            echo = true
            shell = ['/bin/bash','-ue']
            validExitCodes = [0]
            inputs = new LinkedHashMap()
            outputs = new LinkedHashMap()
        }

        configProperties.errorStrategy = ErrorStrategy.TERMINATE
    }


    def methodMissing(String name, Object args) {

        if( args instanceof Object[] ) {

            if( args.size()==1 ) {
                configProperties[ name ] = args[0]
            }
            else {
                configProperties[ name ] = args.toList()
            }

        }
        else {
            configProperties[ name ] = args
        }

        return this
    }

    def void setProperty( String name, Object value ) {
        assert name
        configProperties[name] = value
    }

    def Object propertyMissing(String name) {
        if( name in configProperties ) {
            return configProperties[ name ]
        }

        throw new MissingPropertyException(name,TaskConfig)
    }


    /**
     * Type shortcut to {@code #configProperties.inputs}
     */
    Map<String,DataflowReadChannel> getInputs() {
        return configProperties.inputs
    }

    /**
     * Type shortcut to {@code #configProperties.outputs}
     */
    Map<String,DataflowWriteChannel> getOutputs() {
        return configProperties.outputs
    }

    boolean getEcho() {
        configProperties.echo
    }

    /**
     * Defines one, or more, input channel
     *
     * @param args
     * @return
     */
    TaskConfig input(Map<String,?> args) {
        // wrap by a linked map o guarantee the insertion order

        args?.each { name, value ->
            if ( value instanceof DataflowBroadcast )  {
                inputs.put( name, value.createReadChannel() )
            }
            else if( value instanceof DataflowReadChannel ) {
                inputs.put( name, value )
            }
            // wrap any collections with a DataflowQueue
            else if( value instanceof Collection ) {
                inputs.put( name, Nextflow.channel(value) )
            }
            // wrap any array with a DataflowQueue
            else if ( value && value.class.isArray() ) {
                inputs.put( name, Nextflow.channel(value as List) )
            }
            // wrap a single value with a DataflowVariable
            else {
                inputs.put( name, Nextflow.val(value) )
            }
        }

        return this
    }

    TaskConfig output(String... files) {
        if ( files ) {
            files.each { name -> outputs.put( name, Nextflow.channel() ) }
        }

        return this
    }

    TaskConfig output(Map<String,DataflowWriteChannel> outputs) {
        if ( outputs ) {
            this.outputs.putAll(outputs)
        }

        return this
    }


    TaskConfig errorStrategy( Object value ) {
        assert value

        if( value instanceof ErrorStrategy ) {
            configProperties.errorStrategy = value
        }
        else {
            configProperties.errorStrategy = ErrorStrategy.valueOf(value.toString().toUpperCase())
        }

        return this
    }

    ErrorStrategy getErrorStrategy() {
        configProperties.errorStrategy
    }


    /**
     * The max memory allow to be used to the job
     *
     * @param mem0 The maximum amount of memory expressed as string value,
     *              accepted units are 'B', 'K', 'M', 'G', 'T', 'P'. So for example
     *              {@code maxMemory '100M'}, {@code maxMemory '2G'}, etc.
     */
    TaskConfig maxMemory( Object mem0 ) {
        assert mem0

        configProperties.maxMemory = (mem0 instanceof MemoryUnit) ? mem0 : new MemoryUnit(mem0.toString())

        return this
    }

    /**
     * The max duration time allowed for the job to be executed.
     *
     * @param duration0 The max allowed time expressed as duration string, Accepted units are 'min', 'hour', 'day'.
     *                  For example {@code maxDuration '30 min'}, {@code maxDuration '10 hour'}, {@code maxDuration '2 day'}
     */
    TaskConfig maxDuration( Object duration0 ) {
        assert duration0

        configProperties.maxDuration = (duration0 instanceof Duration) ? duration0 : new Duration(duration0.toString())

        return this
    }

    TaskConfig validExitCodes( Object values ) {

        if( values instanceof List ) {
            configProperties.validExitCodes = values
        }
        else {
            configProperties.validExitCodes = [values]
        }

        return this
    }

    List<Integer> getValidExitCodes() { configProperties.validExitCodes }


}
