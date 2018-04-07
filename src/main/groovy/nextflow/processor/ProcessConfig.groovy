/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import static nextflow.util.CacheHelper.*

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalConfigException
import nextflow.exception.IllegalDirectiveException
import nextflow.executor.BashWrapperBuilder
import nextflow.script.BaseScript
import nextflow.script.DefaultInParam
import nextflow.script.DefaultOutParam
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
/**
 * Holds the process configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessConfig implements Map<String,Object> {

    static final public transient BOOL_YES = ['true','yes','on']

    static final public transient BOOL_NO = ['false','no','off']

    static final public transient LABEL_REGEXP = ~/[a-zA-Z]([a-zA-Z0-9_]*[a-zA-Z0-9]+)?/

    static final public List<String> DIRECTIVES = [
            'afterScript',
            'beforeScript',
            'echo',
            'cache',
            'cpus',
            'container',
            'containerOptions',
            'cleanup',
            'clusterOptions',
            'disk',
            'echo',
            'errorStrategy',
            'executor',
            'ext',
            'instanceType',
            'queue',
            'label',
            'maxErrors',
            'maxForks',
            'maxRetries',
            'memory',
            'module',
            'penv',
            'publishDir',
            'scratch',
            'shell',
            'storeDir',
            'tag',
            'time',
            'validExitStatus',
            // input-output qualifiers
            'file',
            'set',
            'val',
            'each',
            'env',
            'stdin',
            'stdout',
            'stageInMode',
            'stageOutMode'
    ]

    static final List<String> repeatableDirectives = ['label','module','publishDir']

    static final Map<String,Object> DEFAULT_CONFIG = [
            echo: false,
            cacheable: true,
            shell: BashWrapperBuilder.BASH,
            validExitStatus: [0],
            maxRetries: 0,
            maxErrors: -1,
            errorStrategy: ErrorStrategy.TERMINATE
    ]

    @Delegate
    protected final Map<String,Object> configProperties

    private BaseScript ownerScript

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
    ProcessConfig( BaseScript script ) {
        ownerScript = script
        configProperties = new LinkedHashMap()
        configProperties.putAll( DEFAULT_CONFIG )
    }

    /* Only for testing purpose */
    protected ProcessConfig( Map delegate ) {
        configProperties = delegate
    }

    @PackageScope
    static boolean toBool( value )  {
        if( value instanceof Boolean ) {
            return value.booleanValue()
        }

        return value != null && value.toString().toLowerCase() in BOOL_YES
    }

    /**
     * Enable special behavior to allow the configuration object
     * invoking directive method from the process DSL
     *
     * @param value {@code true} enable capture mode, {@code false} otherwise
     * @return The object itself
     */
    @PackageScope
    ProcessConfig enterCaptureMode( boolean value ) {
        this.throwExceptionOnMissingProperty = value
        return this
    }

    private void checkName(String name) {
        if( DIRECTIVES.contains(name) ) return
        if( name == 'when' ) return

        String message = "Unknown process directive: `$name`"
        def alternatives = DIRECTIVES.closest(name)
        if( alternatives.size()==1 ) {
            message += '\n\nDid you mean of these?'
            alternatives.each {
                message += "\n        $it"
            }
        }
        throw new IllegalDirectiveException(message)
    }

    Object invokeMethod(String name, Object args) {
        /*
         * This is need to patch #497 -- what is happening is that when in the config file
         * is defined a directive like `memory`, `cpus`, etc in by using a closure,
         * this closure is interpreted as method definition and it get invoked if a
         * directive with the same name is defined in the process definition.
         * To avoid that the offending property is removed from the map before the method
         * is evaluated.
         */
        if( configProperties.get(name) instanceof Closure )
            configProperties.remove(name)

        this.metaClass.invokeMethod(this,name,args)
    }

    def methodMissing( String name, def args ) {
        checkName(name)

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

    Object put( String name, Object value ) {

        if( name in repeatableDirectives  ) {
            final result = configProperties.get(name)
            configProperties.remove(name)
            this.metaClass.invokeMethod(this, name, value)
            return result
        }
        else {
            return configProperties.put(name,value)
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
    void fakeInput() {
        new DefaultInParam(this)
    }

    void fakeOutput() {
        new DefaultOutParam(this)
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

    protected boolean isValidLabel(String lbl) {
        def p = lbl.indexOf('=')
        if( p==-1 )
            return LABEL_REGEXP.matcher(lbl).matches()

        def left = lbl.substring(0,p)
        def right = lbl.substring(p+1)
        return LABEL_REGEXP.matcher(left).matches() && LABEL_REGEXP.matcher(right).matches()
    }

    ProcessConfig label(String lbl) {
        if( !lbl ) return this
        // -- check that label has a valid syntax
        if( !isValidLabel(lbl) )
            throw new IllegalConfigException("Not a valid process label: $lbl -- Label must consist of alphanumeric characters or '_', must start with an alphabetic character and must end with an alphanumeric character")

        // -- get the current label, it must be a list
        def allLabels = (List)configProperties.get('label')
        if( !allLabels ) {
            allLabels = new ConfigList()
            configProperties.put('label', allLabels)
        }

        // -- avoid duplicates
        if( !allLabels.contains(lbl) )
            allLabels.add(lbl)
        return this
    }

    List<String> getLabels() {
        (List<String>) configProperties.get('label') ?: Collections.emptyList()
    }

    ProcessConfig module( moduleName ) {
        // when no name is provided, just exit
        if( !moduleName )
            return this

        def result = (List)configProperties.module
        if( result == null ) {
            result = new ConfigList()
            configProperties.put('module', result)
        }

        if( moduleName instanceof List )
            result.addAll(moduleName)
        else
            result.add(moduleName)
        return this
    }


    ProcessConfig errorStrategy( value ) {
        configProperties.put('errorStrategy', value)
        return this
    }

    /**
     * Allow the user to specify publishDir directive as a map eg:
     *
     *     publishDir path:'/some/dir', mode: 'copy'
     *
     * @param params
     *      A map representing the the publishDir setting
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig publishDir(Map params ) {
        if( !params )
            return this

        def dirs = (List)configProperties.get('publishDir')
        if( !dirs ) {
            dirs = new ConfigList()
            configProperties.put('publishDir', dirs)
        }

        dirs.add(params)
        return this
    }

    /**
     * Allow the user to specify publishDir directive with a path and a list of named parameters, eg:
     *
     *     publishDir '/some/dir', mode: 'copy'
     *
     * @param params
     *      A map representing the publishDir properties
     * @param target
     *      The target publishDir path
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig publishDir(Map params, target ) {
        params.put('path', target)
        publishDir( params )
    }

    /**
     * Allow the user to specify the publishDir as a string path, eg:
     *
     *      publishDir '/some/dir'
     *
     * @param target
     *      The target publishDir path
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig publishDir( target ) {
        if( target instanceof List ) {
            for( Object item : target ) { publishDir(item) }
        }
        else if( target instanceof Map ) {
            publishDir( target as Map )
        }
        else {
            publishDir([path: target])
        }
        return this
    }
}
