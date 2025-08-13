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

package nextflow.script

import static nextflow.util.CacheHelper.*

import java.util.regex.Pattern

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.ast.NextflowDSLImpl
import nextflow.exception.ConfigParseException
import nextflow.exception.IllegalConfigException
import nextflow.exception.IllegalDirectiveException
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.ConfigList
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskConfig
import nextflow.script.params.CmdEvalParam
import nextflow.script.params.DefaultInParam
import nextflow.script.params.DefaultOutParam
import nextflow.script.params.EachInParam
import nextflow.script.params.EnvInParam
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.InputsList
import nextflow.script.params.OutParam
import nextflow.script.params.OutputsList
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.TupleInParam
import nextflow.script.params.TupleOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam
import nextflow.util.TestOnly

/**
 * Holds the process configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessConfig implements Map<String,Object>, Cloneable {

    static final public transient LABEL_REGEXP = ~/[a-zA-Z]([a-zA-Z0-9_]*[a-zA-Z0-9]+)?/

    static final public List<String> DIRECTIVES = [
            'accelerator',
            'afterScript',
            'arch',
            'array',
            'beforeScript',
            'cache',
            'cleanup',
            'clusterOptions',
            'conda',
            'container',
            'containerOptions',
            'cpus',
            'debug',
            'disk',
            'echo', // deprecated
            'errorStrategy',
            'executor',
            'ext',
            'fair',
            'label',
            'machineType',
            'maxSubmitAwait',
            'maxErrors',
            'maxForks',
            'maxRetries',
            'memory',
            'module',
            'package',
            'penv',
            'pixi',
            'pod',
            'publishDir',
            'queue',
            'resourceLabels',
            'resourceLimits',
            'scratch',
            'secret',
            'shell',
            'spack',
            'stageInMode',
            'stageOutMode',
            'storeDir',
            'tag',
            'time',
            // input-output qualifiers
            'file',
            'val',
            'each',
            'env',
            'stdin',
            'stdout',
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
            debug: false,
            cacheable: true,
            shell: BashWrapperBuilder.BASH,
            maxRetries: 1,
            maxErrors: -1,
            errorStrategy: ErrorStrategy.TERMINATE
    ]

    /**
     * Map instance that holds the process configuration
     */
    @Delegate
    @PackageScope
    protected Map<String,Object> configProperties

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
    protected ProcessConfig( BaseScript script ) {
        ownerScript = script
        configProperties = new LinkedHashMap()
        configProperties.putAll( DEFAULT_CONFIG )
    }

    ProcessConfig( BaseScript script, String name ) {
        this(script)
        this.processName = name
    }

    @TestOnly
    protected ProcessConfig( Map delegate ) {
        configProperties = delegate
    }

    @Override
    ProcessConfig clone() {
        def copy = (ProcessConfig)super.clone()
        copy.@configProperties = new LinkedHashMap<>(configProperties)
        copy.@inputs = inputs.clone()
        copy.@outputs = outputs.clone()
        return copy
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
        if( DIRECTIVES.contains(name) )
            return
        if( name == NextflowDSLImpl.PROCESS_WHEN )
            return
        if( name == NextflowDSLImpl.PROCESS_STUB )
            return

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

    @Override
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

    TaskConfig createTaskConfig() {
        return new TaskConfig(configProperties)
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
     *      A map object modelling the setting defined defined by the user in the nextflow configuration file
     * @param labels
     *      All the labels representing the object holding the configuration setting to apply
     */
    protected void applyConfigSelectorWithLabels(Map<String,?> configDirectives, List<String> labels ) {
        final prefix = 'withLabel:'
        for( String rule : configDirectives.keySet() ) {
            if( !rule.startsWith(prefix) )
                continue
            final pattern = rule.substring(prefix.size()).trim()
            if( !matchesLabels(labels, pattern) )
                continue

            log.debug "Config settings `$rule` matches labels `${labels.join(',')}` for process with name $processName"
            def settings = configDirectives.get(rule)
            if( settings instanceof Map ) {
                applyConfigSettings(settings)
            }
            else if( settings != null ) {
                throw new ConfigParseException("Unknown config settings for process labeled ${labels.join(',')} -- settings=$settings ")
            }
        }
    }

    static boolean matchesLabels( List<String> labels, String pattern ) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()

        final regex = Pattern.compile(pattern)
        for (label in labels) {
            if (regex.matcher(label).matches()) {
                return !isNegated
            }
        }

        return isNegated
    }

    protected void applyConfigSelectorWithName(Map<String,?> configDirectives, String target ) {
        final prefix = 'withName:'
        for( String rule : configDirectives.keySet() ) {
            if( !rule.startsWith(prefix) )
                continue
            final pattern = rule.substring(prefix.size()).trim()
            if( !matchesSelector(target, pattern) )
                continue

            log.debug "Config settings `$rule` matches process $processName"
            def settings = configDirectives.get(rule)
            if( settings instanceof Map ) {
                applyConfigSettings(settings)
            }
            else if( settings != null ) {
                throw new ConfigParseException("Unknown config settings for process with name: $target  -- settings=$settings ")
            }
        }
    }

    static boolean matchesSelector( String name, String pattern ) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()
        return Pattern.compile(pattern).matcher(name).matches() ^ isNegated
    }

    /**
     * Apply the process configuration provided in the nextflow configuration file
     * to the process instance
     *
     * @param configProcessScope The process configuration settings specified
     *      in the configuration file as {@link Map} object
     * @param simpleName The process name
     */
    void applyConfig(Map configProcessScope, String baseName, String simpleName, String fullyQualifiedName) {
        // -- Apply the directives defined in the config object using the`withLabel:` syntax
        final processLabels = this.getLabels() ?: ['']
        this.applyConfigSelectorWithLabels(configProcessScope, processLabels)

        // -- apply setting defined in the config file using the process base name
        this.applyConfigSelectorWithName(configProcessScope, baseName)

        // -- apply setting defined in the config file using the process simple name
        if( simpleName && simpleName!=baseName )
            this.applyConfigSelectorWithName(configProcessScope, simpleName)

        // -- apply setting defined in the config file using the process qualified name (ie. with the execution scope)
        if( fullyQualifiedName && (fullyQualifiedName!=simpleName || fullyQualifiedName!=baseName) )
            this.applyConfigSelectorWithName(configProcessScope, fullyQualifiedName)

        // -- Apply defaults
        this.applyConfigDefaults(configProcessScope)

        // -- check for conflicting settings
        if( this.scratch && this.stageInMode == 'rellink' ) {
            log.warn("Directives `scratch` and `stageInMode=rellink` conflict with each other -- Enforcing default stageInMode for process `$simpleName`")
            this.remove('stageInMode')
        }
    }

    void applyConfigLegacy(Map configProcessScope, String processName) {
        applyConfig(configProcessScope, processName, null, null)
    }


    /**
     * Apply the settings defined in the configuration file to the actual process configuration object
     *
     * @param settings
     *      A map object modelling the setting defined defined by the user in the nextflow configuration file
     */
    protected void applyConfigSettings(Map<String,?> settings) {
        if( !settings )
            return

        for( Entry<String,?> entry: settings ) {
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
    protected void applyConfigDefaults( Map processDefaults ) {
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

    /**
     * Implements the process {@code debug} directive.
     */
    ProcessConfig debug( value ) {
        configProperties.debug = value
        return this
    }

    /**
     * Implements the process {@code echo} directive for backwards compatibility.
     *
     * note: without this method definition {@link BaseScript#echo} will be invoked
     */
    ProcessConfig echo( value ) {
        log.warn1('The `echo` directive has been deprecated - use to `debug` instead')
        configProperties.debug = value
        return this
    }

    /// input parameters

    InParam _in_val( obj ) {
        new ValueInParam(this).bind(obj)
    }

    InParam _in_file( obj ) {
        new FileInParam(this).bind(obj)
    }

    InParam _in_path( Map opts=null, obj ) {
        new FileInParam(this)
                .setPathQualifier(true)
                .setOptions(opts)
                .bind(obj)
    }

    InParam _in_each( obj ) {
        new EachInParam(this).bind(obj)
    }

    InParam _in_tuple( Object... obj ) {
        new TupleInParam(this).bind(obj)
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

    OutParam _out_val( Map opts, Object obj ) {
        new ValueOutParam(this)
                .setOptions(opts)
                .bind(obj)
    }

    OutParam _out_env( Object obj ) {
        new EnvOutParam(this).bind(obj)
    }

    OutParam _out_env( Map opts, Object obj ) {
        new EnvOutParam(this)
                .setOptions(opts)
                .bind(obj)
    }

    OutParam _out_eval(Object obj ) {
        new CmdEvalParam(this).bind(obj)
    }

    OutParam _out_eval(Map opts, Object obj ) {
        new CmdEvalParam(this)
            .setOptions(opts)
            .bind(obj)
    }

    OutParam _out_file( Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' )
            new StdOutParam(this).bind(obj)

        else
            new FileOutParam(this).bind(obj)
    }

    OutParam _out_path( Map opts=null, Object obj ) {
        // note: check that is a String type to avoid to force
        // the evaluation of GString object to a string
        if( obj instanceof String && obj == '-' ) {
            new StdOutParam(this)
                    .setOptions(opts)
                    .bind(obj)
        }
        else {
            new FileOutParam(this)
                    .setPathQualifier(true)
                    .setOptions(opts)
                    .bind(obj)
        }
    }

    OutParam _out_tuple( Object... obj ) {
        new TupleOutParam(this) .bind(obj)
    }

    OutParam _out_tuple( Map opts, Object... obj ) {
        new TupleOutParam(this)
                .setOptions(opts)
                .bind(obj)
    }

    OutParam _out_stdout( Map opts ) {
        new StdOutParam(this)
                .setOptions(opts)
                .bind('-')
    }

    OutParam _out_stdout( obj = null ) {
        def result = new StdOutParam(this).bind('-')
        if( obj ) {
            result.into(obj)
        }
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
        HashMode.of(configProperties.cache) ?: HashMode.DEFAULT()
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

    /**
     * Implements the process {@code label} directive.
     *
     * Note this directive  can be specified (invoked) more than one time in
     * the process context.
     *
     * @param map
     *      The map to be attached to the process.
     * @return
     *      The {@link ProcessConfig} instance itself.
     */
    ProcessConfig resourceLabels(Map<String, Object> map) {
        if( !map )
            return this

        // -- get the current sticker, it must be a Map
        def allLabels = (Map)configProperties.get('resourceLabels')
        if( !allLabels ) {
            allLabels = [:]
        }
        // -- merge duplicates
        allLabels += map
        configProperties.put('resourceLabels', allLabels)
        return this
    }

    Map<String,Object> getResourceLabels() {
        (configProperties.get('resourceLabels') ?: Collections.emptyMap()) as Map<String, Object>
    }

    List<String> getLabels() {
        (List<String>) configProperties.get('label') ?: Collections.<String>emptyList()
    }

    boolean getFair() {
        final value = configProperties.get('fair')
        if( value == null )
            return false
        if( value instanceof Boolean )
            return value

        if( value instanceof Closure )
            throw new IllegalArgumentException("Process directive `fair` cannot be declared in a dynamic manner with a closure")
        else
            throw new IllegalArgumentException("Unexpected value for directive `fair` -- offending value: $value")
    }

    ProcessConfig secret(String name) {
        if( !name )
            return this

        // -- get the current label, it must be a list
        def allSecrets = (List)configProperties.get('secret')
        if( !allSecrets ) {
            allSecrets = new ConfigList()
            configProperties.put('secret', allSecrets)
        }

        // -- avoid duplicates
        if( !allSecrets.contains(name) )
            allSecrets.add(name)
        return this
    }

    List<String> getSecret() {
        (List<String>) configProperties.get('secret') ?: Collections.<String>emptyList()
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
        if( strategy instanceof CharSequence && !ErrorStrategy.isValid(strategy) ) {
            throw new IllegalArgumentException("Unknown error strategy '${strategy}' ― Available strategies are: ${ErrorStrategy.values().join(',').toLowerCase()}")
        }

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

    ProcessConfig accelerator( Map params, value )  {
        if( value instanceof Number ) {
            if( params.limit==null )
                params.limit=value
            else if( params.request==null )
                params.request=value
        }
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
        accelerator(params)
        return this
    }

    ProcessConfig accelerator( value ) {
        if( value instanceof Number )
            configProperties.put('accelerator', [limit: value])
        else if( value instanceof Map )
            configProperties.put('accelerator', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
        return this
    }

    /**
     * Allow user to specify `disk` directive as a value with a list of options, eg:
     *
     *     disk 375.GB, type: 'local-ssd'
     *
     * @param opts
     *      A map representing the disk options
     * @param value
     *      The default disk value
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig disk( Map opts, value )  {
        opts.request = value
        return disk(opts)
    }

    /**
     * Allow user to specify `disk` directive as a value or a list of options, eg:
     *
     *     disk 100.GB
     *     disk request: 375.GB, type: 'local-ssd'
     *
     * @param value
     *      The default disk value or map of options
     * @return
     *      The {@link ProcessConfig} instance itself
     */
    ProcessConfig disk( value ) {
        if( value instanceof Map || value instanceof Closure )
            configProperties.put('disk', value)
        else
            configProperties.put('disk', [request: value])
        return this
    }

    ProcessConfig arch( Map params, value )  {
        if( value instanceof String ) {
            if( params.name==null )
                params.name=value
        }
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
        arch(params)
        return this
    }

    ProcessConfig arch( value ) {
        if( value instanceof String )
            configProperties.put('arch', [name: value])
        else if( value instanceof Map )
            configProperties.put('arch', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
        return this
    }

    int getArray() {
        final value = configProperties.get('array')
        if( value == null )
            return 0
        if( value instanceof Closure )
            throw new IllegalArgumentException("Process directive `array` cannot be declared in a dynamic manner with a closure")
        try {
            final result = value as Integer
            if( result < 0 )
                throw new IllegalArgumentException("Process directive `array` cannot be a negative number")
            if( result == 1 )
                throw new IllegalArgumentException("Process directive `array` should be greater than 1")
            return result
        }
        catch( NumberFormatException e ) {
            throw new IllegalArgumentException("Process directive `array` should be an integer greater than 1 -- offending value: '$value'", e)
        }
    }

    private static final List<String> VALID_RESOURCE_LIMITS = List.of('cpus', 'memory', 'disk', 'time')

    ProcessConfig resourceLimits( Map entries ) {
        for( entry in entries )
            if( entry.key !in VALID_RESOURCE_LIMITS )
                throw new IllegalArgumentException("Not a valid directive in `resourceLimits`: $entry.key")

        configProperties.put('resourceLimits', entries)
        return this
    }

}
