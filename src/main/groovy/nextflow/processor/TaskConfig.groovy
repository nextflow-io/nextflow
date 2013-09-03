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
import groovy.util.logging.Slf4j
import nextflow.script.BaseScript
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Holds the task configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskConfig implements Map {

    @Delegate
    protected final Map configProperties

    private final BaseScript ownerScript

    /**
     * Initialize the taskConfig object with the defaults values
     */
    TaskConfig( BaseScript script ) {

        ownerScript = script

        configProperties = new LinkedHashMap()
        configProperties.with {
            echo = false
            cacheable = true
            shell = ['/bin/bash','-ue']
            validExitCodes = [0]
            inputs = new InputsList()
            outputs = new OutputsList()
        }

        configProperties.errorStrategy = ErrorStrategy.TERMINATE
    }

    protected TaskConfig( Map delegate ) {
        configProperties = delegate
    }

    protected TaskConfig( TaskConfig cfg ) {
        configProperties = cfg
        ownerScript = cfg.@ownerScript
        log.debug "TaskConfig >> ownerScript: $ownerScript"
    }

    def boolean containsKey(String name) {
        return configProperties.containsKey(name)
    }


    def methodMissing( String name, def args) {

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

    /**
     * Type shortcut to {@code #configProperties.inputs}
     */
    InputsList getInputs() {
        return configProperties.inputs
    }

    /**
     * Type shortcut to {@code #configProperties.outputs}
     */
    OutputsList getOutputs() {
        return configProperties.outputs
    }

    boolean getEcho() {
        configProperties.echo
    }

    void setEcho( Object value ) {
        if( value instanceof Boolean ) {
            configProperties.echo = value.booleanValue()
        }
        else if( value ?.toString()?.toLowerCase() in ['true','yes','on']) {
            configProperties.echo = true
        }
        else {
            configProperties.echo = false
        }
    }

    TaskConfig echo( def value ) {
        setEcho(value)
        return this
    }


    /**
     * Defines a task 'input'
     *
     * @param attributes
     * @return
     */
    def TaskConfig input( Map attributes )  {

        configProperties.inputs << InParam.create(attributes)

        return this
    }

    /**
     * Defines a special *dummy* input parameter, when no inputs are
     * provided by the user for the current task
     */
    def void noInput() {
        input( val:'$', from: true )
    }

    /**
     * Defines a task 'output'. An output is defined through a map holding its
     * properties, using the following pattern:
     * <p>
     * <pre>
     *     file: <file name>, into: <receiving channel>
     * </pre>
     *
     *
     * @param attributes
     * @return
     */
    def TaskConfig output( Map attributes ) {

        // add it to the list out outputs
        configProperties.outputs << OutParam.create(attributes, ownerScript)

        return this
    }

    def TaskConfig stdin(channel) {
        input( file:'-', from:channel )
    }

    def TaskConfig stdout(channel) {
        output( file:'-', into:channel )
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

    Duration getMaxDuration() { configProperties.maxDuration }

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
