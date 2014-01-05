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
import nextflow.ast.ProcessVarRef
import nextflow.script.BaseScript
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.InputsList
import nextflow.script.OutputsList
import nextflow.script.SharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.script.ValueSharedParam
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
            undef = false
            cacheable = true
            shell = ['/bin/bash','-ue']
            validExitStatus = [0]
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
        log.trace "TaskConfig >> ownerScript: $ownerScript"
    }

    private boolean toBool( value )  {
        if( value instanceof Boolean ) {
            configProperties.echo = value.booleanValue()
        }
        else if( value != null && value.toString().toLowerCase() in ['true','yes','on'] ) {
            configProperties.echo = true
        }
        else {
            configProperties.echo = false
        }
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

    List getSharedDefs () {
        configProperties.inputs.findAll { it instanceof SharedParam }
    }

    boolean getEcho() {
        configProperties.echo
    }

    boolean getMerge() {
        configProperties.merge?.toString()?.toLowerCase() in ['true','yes','on']
    }

    void setEcho( Object value ) {
        configProperties.echo = toBool(value)
    }

    TaskConfig echo( def value ) {
        setEcho(value)
        return this
    }

    void setUndef( def value ) {
        configProperties.undef = toBool(value)
    }

    TaskConfig undef( def value ) {
        configProperties.undef = toBool(value)
        return this
    }

    boolean getUndef() {
        configProperties.undef
    }

    /**
     * Define a *file* parameter which can be used both input or output declaration
     *
     * @param obj The file name
     * @param direction The parameter direction, i.e. {@code in} or {@code out}
     * @return
     *
     */
    def FileInParam _in_file( boolean isVariable, Object obj ) {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new FileInParam(ownerScript, param)
        configProperties.inputs << result

        return result
    }

    def ValueInParam _in_val( boolean isVariable, Object obj ) {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new ValueInParam(ownerScript, param)
        configProperties.inputs << result

        return result
    }

    EnvInParam _in_env( boolean isVariable, Object obj ) {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new EnvInParam(ownerScript, param)
        configProperties.inputs << result

        result
    }

    EachInParam _in_each( boolean isVariable, Object obj ) {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new EachInParam(ownerScript,param)
        configProperties.inputs << result

        return result
    }

    def FileOutParam _out_file( String name ) {

        def result = name == '-' ? new StdOutParam(ownerScript) : new FileOutParam(ownerScript,name)
        configProperties.outputs << result

        result
    }

    def ValueOutParam _out_val( String name ) {

        def result = new ValueOutParam(ownerScript,name)
        configProperties.outputs << result

        result
    }


    def ValueSharedParam _share_val( boolean isVariable, def obj )  {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new ValueSharedParam(ownerScript, param)
        configProperties.inputs << result

        return result
    }

    def FileSharedParam _share_file( boolean isVariable, def obj )  {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new FileSharedParam(ownerScript, param)
        configProperties.inputs << result

        return result
    }



    StdInParam stdin( boolean isVariable, def obj ) {

        def param = isVariable ? new ProcessVarRef(obj.toString()) : obj
        def result = new StdInParam(ownerScript, param)

        configProperties.inputs << result

        result
    }

    StdOutParam stdout( def channelRef = null ) {

        def result = new StdOutParam(ownerScript)
        if( channelRef != null )
            result.to(channelRef)

        configProperties.outputs << result

        result
    }

    /**
     * Defines a special *dummy* input parameter, when no inputs are
     * provided by the user for the current task
     */
    def void noInput() {
        _in_val(false,'$') ._as(true)
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

    TaskConfig validExitStatus( Object values ) {

        if( values instanceof List ) {
            configProperties.validExitStatus = values
        }
        else {
            configProperties.validExitStatus = [values]
        }

        return this
    }

    List<Integer> getValidExitStatus() { configProperties.validExitStatus }



}
