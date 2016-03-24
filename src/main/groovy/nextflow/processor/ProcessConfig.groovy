/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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
import static nextflow.util.CacheHelper.HashMode

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.executor.BashWrapperBuilder
import nextflow.script.BaseScript
import nextflow.script.EachInParam
import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.InParam
import nextflow.script.InputsList
import nextflow.script.OutParam
import nextflow.script.OutputsList
import nextflow.script.SetInParam
import nextflow.script.SetOutParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import nextflow.util.ReadOnlyMap
/**
 * Holds the process configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessConfig implements Map<String,Object> {

    static final transient BOOL_YES = ['true','yes','on']

    static final transient BOOL_NO = ['false','no','off']

    @Delegate
    protected final Map<String,Object> configProperties

    private final BaseScript ownerScript

    private boolean throwExceptionOnMissingProperty

    private inputs = new InputsList()

    private outputs = new OutputsList()

    /**
     * Initialize the taskConfig object with the defaults values
     *
     * @param script The owner {@code BaseScript} configuration object
     * @param important The values specified by this map won't be overridden by attributes
     *      having the same name defined at the task level
     */
    ProcessConfig( BaseScript script, Map important = null ) {

        ownerScript = script

        // parse the attribute as List before adding it to the read-only list
        if( important?.containsKey('module') ) {
            important.module = parseModule(important.module)
        }

        configProperties = important ? new ReadOnlyMap(important) : new LinkedHashMap()
        configProperties.echo = false
        configProperties.cacheable = true
        configProperties.shell = BashWrapperBuilder.BASH
        configProperties.validExitStatus = [0]
        configProperties.maxRetries = 1
        configProperties.maxErrors = 3
        configProperties.errorStrategy = ErrorStrategy.TERMINATE
    }

    /* Only for testing purpose */
    protected ProcessConfig( Map delegate ) {
        configProperties = delegate
    }

    protected ProcessConfig( ProcessConfig cfg ) {
        configProperties = cfg
        ownerScript = cfg.@ownerScript
        log.trace "TaskConfig >> ownerScript: $ownerScript"
    }

    @PackageScope
    static boolean toBool( value )  {
        if( value instanceof Boolean ) {
            return value.booleanValue()
        }

        return value != null && value.toString().toLowerCase() in BOOL_YES
    }

    @PackageScope
    ProcessConfig throwExceptionOnMissingProperty( boolean value ) {
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
            case 'inputs':
                return getInputs()

            case 'outputs':
                return getOutputs()

            case 'cacheable':
                return isCacheable()

            case 'ext':
                if( !configProperties.containsKey('ext') ) {
                    configProperties.put('ext', new HashMap())
                }
                return configProperties.get('ext')

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

    @PackageScope
    TaskConfig createTaskConfig() {
        new TaskConfig(configProperties)
    }

    /**
     * Type shortcut to {@code #configProperties.inputs}
     */
    InputsList getInputs() {
        inputs
    }

    /**
     * Type shortcut to {@code #configProperties.outputs}
     */
    OutputsList getOutputs() {
        outputs
    }


    /*
     * note: without this method definition {@link BaseScript#echo} will be invoked
     */
    ProcessConfig echo( value ) {
        configProperties.echo = value
        return this
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
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(this).bind(obj)

        else
            new FileOutParam(this).bind(obj)
    }

    OutParam _out_set( Object... obj ) {
        new SetOutParam(this) .bind(obj)
    }

    OutParam _out_stdout( obj = null ) {
        def result = new StdOutParam(this).bind('-')
        if( obj ) result.into(obj)
        result
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

    ProcessConfig module( moduleName ) {
        // when no name is provided, just exit
        if( !moduleName )
            return this

        def list = parseModule(moduleName, configProperties.module)
        configProperties.put('module', list)
        return this
    }

    @PackageScope
    static List parseModule( value, current = null) {

        // if no modules list exist create it
        List copy

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
            copy.add( value )

        return copy
    }

}
