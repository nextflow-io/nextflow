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
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.executor.BashWrapperBuilder
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
class TaskConfig implements Map<String,Object> {

    static final transient BOOL_YES = ['true','yes','on']

    static final transient BOOL_NO = ['false','no','off']

    @Delegate
    protected final Map<String,Object> configProperties

    private final BaseScript ownerScript

    private boolean throwExceptionOnMissingProperty

    /**
     * Initialize the taskConfig object with the defaults values
     *
     * @param script The owner {@code BaseScript} configuration object
     * @param important The values specified by this map won't be overridden by attributes
     *      having the same name defined at the task level
     */
    TaskConfig( BaseScript script, Map important = null ) {

        ownerScript = script

        // parse the attribute as List before adding it to the read-only list
        if( important?.containsKey('module') ) {
            important.module = parseModuleString(important.module)
        }

        configProperties = important ? new ReadOnlyMap(important) : new LinkedHashMap()
        configProperties.with {
            echo = false
            undef = false
            cacheable = true
            shell = (BashWrapperBuilder.BASH)
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
            return value.booleanValue()
        }

        return value != null && value.toString().toLowerCase() in BOOL_YES
    }

    @PackageScope
    TaskConfig throwExceptionOnMissingProperty( boolean value ) {
        this.throwExceptionOnMissingProperty = value
        return this
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

    def getProperty( String name ) {

        switch( name ) {
            case 'cacheable':
                return isCacheable()

            case 'errorStrategy':
                return getErrorStrategy()

            case 'shell':
                return getShell()

            case 'time':
                return getTime()

            case 'memory':
                return getMemory()

            default:
                if( configProperties.containsKey(name) )
                    return configProperties.get(name)
                else if( throwExceptionOnMissingProperty )
                    throw new MissingPropertyException("Unknown variable '$name'", name, null)
                else
                    return null
        }

    }

    @PackageScope
    BaseScript getOwnerScript() { ownerScript }

    /**
     * Type shortcut to {@code #configProperties.inputs}
     */
    InputsList getInputs() {
        configProperties.inputs
    }

    /**
     * Type shortcut to {@code #configProperties.outputs}
     */
    OutputsList getOutputs() {
        configProperties.outputs
    }

    List getSharedDefs () {
        configProperties.inputs.findAll { it instanceof SharedParam }
    }

    boolean getEcho() {
        configProperties.echo
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
     * @param value The maximum amount of memory expressed as string value,
     *              accepted units are 'B', 'K', 'M', 'G', 'T', 'P'. So for example
     *              {@code maxMemory '100M'}, {@code maxMemory '2G'}, etc.
     */
    @Deprecated
    TaskConfig maxMemory( Object value ) {
        log.warn("Directive 'maxMemory' has been deprecated. Use 'memory' instead.")
        configProperties.memory = value
        return this
    }

    /**
     * The max duration time allowed for the job to be executed.
     *
     * @param value The max allowed time expressed as duration string, Accepted units are 'min', 'hour', 'day'.
     *                  For example {@code maxDuration '30 min'}, {@code maxDuration '10 hour'}, {@code maxDuration '2 day'}
     */
    @Deprecated
    TaskConfig maxDuration( Object value ) {
        log.warn("Directive 'maxDuration' has been deprecated. Use 'time' instead.")
        configProperties.time = value
        return this
    }


    MemoryUnit getMemory() {
        def value = configProperties.memory

        if( !value )
            return null

        if( value instanceof MemoryUnit )
            return value

        try {
            new MemoryUnit(value.toString().trim())
        }
        catch( Exception e ) {
            throw new AbortOperationException("Not a valid 'memory' value in process definition: $value")
        }
    }

    Duration getTime() {
        def value = configProperties.time

        if( !value )
            return null

        if( value instanceof Duration )
            return value

        try {
            new Duration(value.toString().trim())
        }
        catch( Exception e ) {
            throw new AbortOperationException("Not a valid `time` value in process definition: $value")
        }
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

    List<Integer> getValidExitStatus() {
        (List<Integer>)configProperties.validExitStatus
    }

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

    TaskConfig module( moduleName ) {
        // when no name is provided, just exit
        if( !moduleName )
            return this

        def list = parseModuleString(moduleName, configProperties.module)
        configProperties.put('module', list)
        return this
    }

    private List<String> parseModuleString( value, current = null) {

        // if no modules list exist create it
        List<String> copy

        // normalize the current value to a list
        // note: any value that is not a list is discarded
        if( current instanceof List )
            copy = new ArrayList<>(current)
        else
            copy = []

        // parse the module list
        if( value instanceof List )
            copy.addAll(value)
        else if( value instanceof String && value.contains(':'))
            for( String it : value.split(':') ) { copy.add(it) }
        else
            copy.add( value.toString() )

        return copy
    }

    List<String> getModule() {
        def result = configProperties.module
        if( result instanceof String ) {
            result = configProperties.module = parseModuleString(result)
        }
        (List<String>) result
    }

    List<String> getShell() {
        final value = configProperties.shell
        if( !value )
            return BashWrapperBuilder.BASH

        if( value instanceof List )
            return value

        if( value instanceof CharSequence )
            return [ value.toString() ]

        throw new IllegalArgumentException("Not a valid 'shell' configuration value: ${value}")
    }

}
