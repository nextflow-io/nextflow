/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.processor

import static nextflow.util.CacheHelper.*

import java.util.regex.Pattern

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.exception.ConfigParseException
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

    static final public transient LABEL_REGEXP = ~/[a-zA-Z]([a-zA-Z0-9_]*[a-zA-Z0-9]+)?/

    static final public List<String> DIRECTIVES = [
            'afterScript',
            'beforeScript',
            'echo',
            'cache',
            'conda',
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
            'pod',
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

    /**
     * Names of directives that can be used more than once in the process definition
     */
    @PackageScope
    static final List<String> repeatableDirectives = ['label','module','pod','publishDir']

    /**
     * Default directives values
     */
    @PackageScope
    static final Map<String,Object> DEFAULT_CONFIG = [
            echo: false,
            cacheable: true,
            shell: BashWrapperBuilder.BASH,
            validExitStatus: [0],
            maxRetries: 0,
            maxErrors: -1,
            errorStrategy: ErrorStrategy.TERMINATE
    ]

    /**
     * Map instance that holds the process configuration
     */
    @Delegate
    @PackageScope
    protected final Map<String,Object> configProperties

    /**
     * Reference to the main script object
     */
    private BaseScript ownerScript

    /**
     * Name of the process to which this config object is associated
     */
    private String processName

    /**
     * When {@code true} a {@link MissingPropertyException} is thrown when
     * trying to access a property not existing
     */
    private boolean throwExceptionOnMissingProperty

    /**
     * List of process input definitions
     */
    private inputs = new InputsList()

    /**
     * List of process output definitions
     */
    private outputs = new OutputsList()

    /**
     * Initialize the taskConfig object with the defaults values
     *
     * @param script The owner {@code BaseScript} configuration object
     */
    ProcessConfig( BaseScript script ) {
        ownerScript = script
        configProperties = new LinkedHashMap()
        configProperties.putAll( DEFAULT_CONFIG )
    }

    /* Only for testing purpose */
    @PackageScope
    ProcessConfig( Map delegate ) {
        configProperties = delegate
    }

    /**
     * Define the name of the process to which this config object is associated
     *
     * @param name String representing the name of the process to which this config object is associated
     * @return The {@link ProcessConfig} instance itself
     */
    @PackageScope
    ProcessConfig setProcessName( String name ) {
        this.processName = name
        return this
    }


    /**
     * Enable special behavior to allow the configuration object
     * invoking directive method from the process DSL
     *
     * @param value {@code true} enable capture mode, {@code false} otherwise
     * @return The object itself
     */
    @PackageScope
    ProcessConfig throwExceptionOnMissingProperty( boolean value ) {
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

    Object getProperty( String name ) {

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
     * Apply the settings defined in the configuration file for the given annotation label, for example:
     *
     * ```
     * process {
     *     withLabel: foo {
     *         cpus = 1
     *         memory = 2.gb
     *     }
     * }
     * ```
     *
     * @param configDirectives
     *      A map object modelling the setting defined defined by the user in the nextflow configuration file     *
     * @param category
     *      The type of annotation either {@code withLabel:} or {@code withName:}
     * @param processLabel
     *      A specific label representing the object holding the configuration setting to apply
     */
    @PackageScope
    void applyConfigForLabel( Map<String,?> configDirectives, String category, String processLabel ) {
        assert category in ['withLabel:','withName:']
        assert processName != null
        assert processLabel != null

        for( String rule : configDirectives.keySet() ) {
            if( !rule.startsWith(category) )
                continue
            final isLabel = category=='withLabel:'
            final pattern = rule.substring(category.size()).trim()
            final target = isLabel ? processLabel : processName
            if( !matchesSelector(target, pattern) )
                continue

            log.debug "Config settings `$rule` matches ${isLabel ? "label `$target` for process with name $processName" : "process $processName"}"
            def settings = configDirectives.get(rule)
            if( settings instanceof Map ) {
                applyConfigSettings(settings)
            }
            else if( settings != null ) {
                throw new ConfigParseException("Unknown config settings for label `$processLabel` -- settings=$settings ")
            }
        }
    }

    static boolean matchesSelector( String target, String pattern ) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()
        return Pattern.compile(pattern).matcher(target).matches() ^ isNegated
    }

    /**
     * Apply the settings defined in the configuration file to the actual process configuration object
     *
     * @param settings
     *      A map object modelling the setting defined defined by the user in the nextflow configuration file
     */
    @PackageScope
    void applyConfigSettings(Map<String,?> settings) {
        if( !settings )
            return

        for( Entry<String,?> entry: settings ) {
            if( entry.key.startsWith('$'))
                continue

            if( entry.key.startsWith("withLabel:") || entry.key.startsWith("withName:"))
                continue

            if( !DIRECTIVES.contains(entry.key) )
                log.warn "Unknown directive `$entry.key` for process `$processName`"

            if( entry.key == 'params' ) // <-- patch issue #242
                continue

            if( entry.key == 'ext' ) {
                if( this.getProperty('ext') instanceof Map ) {
                    // update missing 'ext' properties found in 'process' scope
                    def ext = this.getProperty('ext') as Map
                    entry.value.each { String k, v -> ext[k] = v }
                }
                continue
            }

            this.put(entry.key,entry.value)
        }
    }

    /**
     * Apply the process settings defined globally in the process config scope
     *
     * @param processDefaults
     *      A map object representing the setting to be applied to the process
     *      (provided it does not already define a different value for
     *      the same config setting).
     *
     */
    @PackageScope
    void applyConfigDefaults( Map processDefaults ) {
        for( String key : processDefaults.keySet() ) {
            if( key == 'params' )
                continue
            final value = processDefaults.get(key)
            final current = this.getProperty(key)
            if( key == 'ext' ) {
                if( value instanceof Map && current instanceof Map ) {
                    final ext = current as Map
                    value.each { k,v -> if(!ext.containsKey(k)) ext.put(k,v) }
                }
            }
            else if( !this.containsKey(key) || (DEFAULT_CONFIG.containsKey(key) && current==DEFAULT_CONFIG.get(key)) ) {
                this.put(key, value)
            }
        }
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

        if( value instanceof String && value in Const.BOOL_NO )
            return false

        return true
    }

    HashMode getHashMode() {
        HashMode.of(configProperties.cache) ?: HashMode.STANDARD
    }

    protected boolean isValidLabel(String lbl) {
        def p = lbl.indexOf('=')
        if( p==-1 )
            return LABEL_REGEXP.matcher(lbl).matches()

        def left = lbl.substring(0,p)
        def right = lbl.substring(p+1)
        return LABEL_REGEXP.matcher(left).matches() && LABEL_REGEXP.matcher(right).matches()
    }

    /**
     * Implements the process {@code label} directive.
     *
     * Note this directive  can be specified (invoked) more than one time in
     * the process context.
     *
     * @param lbl
     *      The label to be attached to the process.
     * @return
     *      The {@link ProcessConfig} instance itself.
     */
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

    /**
     * Implements the process {@code module} directive.
     *
     * See also http://modules.sourceforge.net
     *
     * @param moduleName
     *      The module name to be used to execute the process.
     * @return
     *      The {@link ProcessConfig} instance itself.
     */
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

    /**
     * Implements the {@code errorStrategy} directive
     *
     * @see ErrorStrategy
     *
     * @param strategy
     *      A string representing the error strategy to be used.
     * @return
     *      The {@link ProcessConfig} instance itself.
     */
    ProcessConfig errorStrategy( strategy ) {
        configProperties.put('errorStrategy', strategy)
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
    ProcessConfig publishDir(Map params) {
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
    ProcessConfig publishDir(Map params, target) {
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

    private static final Map POD_OPTIONS = [secret:String, mountPath:String, envName: String]

    /**
     * Allow use to specify K8s `pod` options
     *
     * @param entry
     *      A map object representing pod config options
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig pod( Map entry ) {

        if( !entry )
            return this

        def allOptions = (List)configProperties.get('pod')
        if( !allOptions ) {
            allOptions = new ConfigList()
            configProperties.put('pod', allOptions)
        }

        allOptions.add(entry)
        return this

    }
}
