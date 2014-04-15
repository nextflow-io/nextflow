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

package nextflow.processor

import groovy.util.logging.Slf4j
import nextflow.script.BaseScript
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.InParam
import nextflow.script.InputsList
import nextflow.script.OutParam
import nextflow.script.OutputsList
import nextflow.script.SetInParam
import nextflow.script.SetOutParam
import nextflow.script.SharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.script.ValueSharedParam
import nextflow.util.Duration
import nextflow.util.HashMode
import nextflow.util.MemoryUnit
import nextflow.util.ReadOnlyMap

/**
 * Holds the task configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class TaskConfig implements Map {

    private static transient BOOL_YES = ['true','yes','on']

    private static transient BOOL_NO = ['false','no','off']

    @Delegate
    protected final Map configProperties

    private final BaseScript ownerScript

    /**
     * Initialize the taskConfig object with the defaults values
     *
     * @param script The owner {@code BaseScript} configuration object
     * @param important The values specified by this map won't be overridden by attributes
     *      having the same name defined at the task level
     */
    TaskConfig( BaseScript script, Map important = null ) {

        ownerScript = script

        configProperties = important ? new ReadOnlyMap(important) : new LinkedHashMap()
        configProperties.with {
            echo = false
            undef = false
            cacheable = true
            shell = ['/bin/bash','-ue']
            validExitStatus = [0]
            inputs = new InputsList()
            outputs = new OutputsList()
            maxRetries = 1
            maxErrors = 3
        }

        configProperties.errorStrategy = ErrorStrategy.TERMINATE
    }

    /* Only for testing purpose */
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
        else if( value != null && value.toString().toLowerCase() in BOOL_YES ) {
            configProperties.echo = true
        }
        else {
            configProperties.echo = false
        }
    }

    def boolean containsKey(String name) {
        return configProperties.containsKey(name)
    }

    def methodMissing( String name, def args ) {

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

    public getProperty( String name ) {
        switch( name ) {
            case 'cacheable':
                return isCacheable()

            case 'errorStrategy':
                return getErrorStrategy()

            default:
                return configProperties.get(name)
        }
    }


    @groovy.transform.PackageScope
    BaseScript getOwnerScript() { ownerScript }

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
        configProperties.merge?.toString() in BOOL_YES
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

    /// input parameters

    InParam _in_val( obj ) {
        new ValueInParam(this).bind(obj)
    }

    InParam _in_file( obj ) {
        new FileInParam(this).bind(obj)
    }

    InParam _in_each( obj ) {
        new EachInParam(this).bind(obj)
    }

    InParam _in_set( Object... obj ) {
        new SetInParam(this).bind(obj)
    }

    InParam _in_stdin( obj = null ) {
        def result = new StdInParam(this)
        if( obj ) result.bind(obj)
        result
    }

    InParam _in_env( obj ) {
        new EnvInParam(this).bind(obj)
    }


    /// output parameters

    OutParam _out_val( Object obj ) {
        new ValueOutParam(this).bind(obj)
    }


    OutParam _out_file( Object obj ) {
        def result = obj == '-' ? new StdOutParam(this) : new FileOutParam(this)
        result.bind(obj)
    }

    OutParam _out_set( Object... obj ) {
        new SetOutParam(this) .bind(obj)
    }

    OutParam _out_stdout( obj = null ) {
        def result = new StdOutParam(this).bind('-')
        if( obj ) result.into(obj)
        result
    }

    /// shared parameters

    SharedParam _share_val( def obj )  {
        new ValueSharedParam(this).bind(obj) as SharedParam
    }

    SharedParam _share_file( def obj )  {
        new FileSharedParam(this).bind(obj) as SharedParam
    }




    /**
     * Defines a special *dummy* input parameter, when no inputs are
     * provided by the user for the current task
     */
    def void fakeInput() {
        new ValueInParam(this).from(true).bind('$')
    }

    def void fakeOutput( target ) {
        new StdOutParam(this).bind('-').into(target)
    }


    ErrorStrategy getErrorStrategy() {
        switch( configProperties.errorStrategy ) {
            case CharSequence:
                return configProperties.errorStrategy.toUpperCase() as ErrorStrategy
            case ErrorStrategy:
                return configProperties.errorStrategy
            case null:
                return null
            default:
                throw new IllegalArgumentException("Not a valid 'ErrorStrategy' value: ${configProperties.errorStrategy}")
        }
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

    boolean isCacheable() {
        def value = configProperties.cache
        if( value == null )
            return true

        if( value instanceof Boolean )
            return value

        if( value instanceof String && value in BOOL_NO )
            return false

        return true
    }

    HashMode getHashMode() {
        configProperties.cache == 'deep' ? HashMode.DEEP : HashMode.STANDARD
    }

    int getMaxRetries() {
        configProperties.maxRetries ? configProperties.maxRetries as int : 0
    }

    int getMaxErrors() {
        configProperties.maxErrors ? configProperties.maxErrors as int : 0
    }


}
