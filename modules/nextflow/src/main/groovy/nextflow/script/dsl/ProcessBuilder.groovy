/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.script.dsl

import java.util.regex.Pattern

import groovy.util.logging.Slf4j
import nextflow.ast.NextflowDSLImpl
import nextflow.exception.IllegalConfigException
import nextflow.exception.IllegalDirectiveException
import nextflow.exception.ScriptRuntimeException
import nextflow.processor.ErrorStrategy
import nextflow.script.ProcessInputs
import nextflow.script.ProcessOutputs
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ProcessDef
import nextflow.util.LazyList

/**
 * Builder for {@link ProcessDef}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class ProcessBuilder {

    static final List<String> DIRECTIVES = [
            'accelerator',
            'afterScript',
            'arch',
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
            'maxErrors',
            'maxForks',
            'maxRetries',
            'maxSubmitAwait',
            'memory',
            'module',
            'penv',
            'pod',
            'publishDir',
            'queue',
            'resourceLabels',
            'scratch',
            'secret',
            'shell',
            'spack',
            'stageInMode',
            'stageOutMode',
            'storeDir',
            'tag',
            'time'
    ]

    protected BaseScript ownerScript

    protected String processName

    protected BodyDef body

    protected ProcessConfig config

    ProcessBuilder(BaseScript ownerScript, String processName) {
        this.ownerScript = ownerScript
        this.processName = processName
        this.config = new ProcessConfig(ownerScript, processName)
    }

    ProcessBuilder(ProcessConfig config) {
        this.ownerScript = config.getOwnerScript()
        this.processName = config.getProcessName()
        this.config = config
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
        if( config.get(name) instanceof Closure )
            config.remove(name)

        this.metaClass.invokeMethod(this,name,args)
    }

    def methodMissing( String name, def args ) {
        checkName(name)

        if( args instanceof Object[] )
            config.put(name, args.size()==1 ? args[0] : args.toList())
        else
            config.put(name, args)
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
            message += '\n\nDid you mean one of these?'
            alternatives.each {
                message += "\n        $it"
            }
        }
        throw new IllegalDirectiveException(message)
    }

    /// DIRECTIVES

    void accelerator( Map params, value )  {
        if( value instanceof Number ) {
            if( params.limit==null )
                params.limit=value
            else if( params.request==null )
                params.request=value
        }
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
        accelerator(params)
    }

    void accelerator( value ) {
        if( value instanceof Number )
            config.put('accelerator', [limit: value])
        else if( value instanceof Map )
            config.put('accelerator', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
    }

    void arch( Map params, value )  {
        if( value instanceof String ) {
            if( params.name==null )
                params.name=value
        }
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
        arch(params)
    }

    void arch( value ) {
        if( value instanceof String )
            config.put('arch', [name: value])
        else if( value instanceof Map )
            config.put('arch', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
    }

    void debug(boolean value) {
        config.debug = value
    }

    /**
     * Implements the {@code disk} directive, e.g.:
     *
     *     disk 375.GB, type: 'local-ssd'
     *
     * @param opts
     * @param value
     */
    void disk( Map opts, value )  {
        opts.request = value
        disk(opts)
    }

    /**
     * Implements the {@code disk} directive, e.g.:
     *
     *     disk 100.GB
     *     disk request: 375.GB, type: 'local-ssd'
     *
     * @param value
     */
    void disk( value ) {
        if( value instanceof Map || value instanceof Closure )
            config.put('disk', value)
        else
            config.put('disk', [request: value])
    }

    /**
     * Implements the {@code echo} directive for backwards compatibility.
     *
     * note: without this method definition {@link BaseScript#echo} will be invoked
     *
     * @param value
     */
    void echo( value ) {
        log.warn1('The `echo` directive has been deprecated - use `debug` instead')
        config.put('debug', value)
    }

    /**
     * Implements the {@code errorStrategy} directive.
     *
     * @param strategy
     */
    void errorStrategy( CharSequence strategy ) {
        if( !ErrorStrategy.isValid(strategy) )
            throw new IllegalArgumentException("Unknown error strategy '${strategy}' ― Available strategies are: ${ErrorStrategy.values().join(',').toLowerCase()}")

        config.put('errorStrategy', strategy)
    }

    /**
     * Implements the {@code label} directive.
     *
     * This directive can be specified (invoked) more than once in
     * the process definition.
     *
     * @param lbl
     */
    void label(String lbl) {
        if( !lbl ) return

        // -- check that label has a valid syntax
        if( !isValidLabel(lbl) )
            throw new IllegalConfigException("Not a valid process label: $lbl -- Label must consist of alphanumeric characters or '_', must start with an alphabetic character and must end with an alphanumeric character")

        // -- get the current label, it must be a list
        def allLabels = (List)config.get('label')
        if( !allLabels ) {
            allLabels = new LazyList()
            config.put('label', allLabels)
        }

        // -- avoid duplicates
        if( !allLabels.contains(lbl) )
            allLabels.add(lbl)
    }

    private static final Pattern LABEL_REGEXP = ~/[a-zA-Z]([a-zA-Z0-9_]*[a-zA-Z0-9]+)?/

    protected static boolean isValidLabel(String lbl) {
        def p = lbl.indexOf('=')
        if( p==-1 )
            return LABEL_REGEXP.matcher(lbl).matches()

        def left = lbl.substring(0,p)
        def right = lbl.substring(p+1)
        return LABEL_REGEXP.matcher(left).matches() && LABEL_REGEXP.matcher(right).matches()
    }

    /**
     * Implements the {@code module} directive.
     *
     * See also http://modules.sourceforge.net
     *
     * @param value
     */
    void module( String value ) {
        if( !value ) return

        def result = (List)config.module
        if( result == null ) {
            result = new LazyList()
            config.put('module', result)
        }

        result.add(value)
    }

    /**
     * Implements the {@code pod} directive.
     *
     * @param entry
     */
    void pod( Map entry ) {
        if( !entry ) return

        def allOptions = (List)config.get('pod')
        if( !allOptions ) {
            allOptions = new LazyList()
            config.put('pod', allOptions)
        }

        allOptions.add(entry)
    }

    /**
     * Implements the {@code publishDir} directive as a map eg:
     *
     *     publishDir path: '/some/dir', mode: 'copy'
     *
     * This directive can be specified (invoked) multiple times in
     * the process definition.
     *
     * @param params
     */
    void publishDir(Map params) {
        if( !params ) return

        def dirs = (List)config.get('publishDir')
        if( !dirs ) {
            dirs = new LazyList()
            config.put('publishDir', dirs)
        }

        dirs.add(params)
    }

    /**
     * Implements the {@code publishDir} directive as a path with named parameters, eg:
     *
     *     publishDir '/some/dir', mode: 'copy'
     *
     * @param params
     * @param path
     */
    void publishDir(Map params, CharSequence path) {
        params.put('path', path)
        publishDir( params )
    }

    /**
     * Implements the {@code publishDir} directive as a string path, eg:
     *
     *     publishDir '/some/dir'
     *
     * @param target
     */
    void publishDir( target ) {
        if( target instanceof List ) {
            for( Object item : target ) { publishDir(item) }
        }
        else if( target instanceof Map ) {
            publishDir( target as Map )
        }
        else {
            publishDir([path: target])
        }
    }

    /**
     * Implements the {@code resourceLabels} directive.
     *
     * This directive can be specified (invoked) multiple times in
     * the process definition.
     *
     * @param map
     */
    void resourceLabels(Map<String, Object> map) {
        if( !map ) return

        // -- get the current sticker, it must be a Map
        def allLabels = (Map)config.get('resourceLabels')
        if( !allLabels ) {
            allLabels = [:]
        }
        // -- merge duplicates
        allLabels += map
        config.put('resourceLabels', allLabels)
    }

    /**
     * Implements the {@code secret} directive.
     *
     * This directive can be specified (invoked) multiple times in
     * the process definition.
     *
     * @param name
     */
    void secret(String name) {
        if( !name ) return

        // -- get the current label, it must be a list
        def allSecrets = (List)config.get('secret')
        if( !allSecrets ) {
            allSecrets = new LazyList()
            config.put('secret', allSecrets)
        }

        // -- avoid duplicates
        if( !allSecrets.contains(name) )
            allSecrets.add(name)
    }

    /// SCRIPT

    ProcessBuilder withInputs(ProcessInputs inputs) {
        config.inputs = inputs
        return this
    }

    ProcessBuilder withOutputs(ProcessOutputs outputs) {
        config.outputs = outputs
        return this
    }

    ProcessBuilder withBody(Closure closure, String section, String source='', List values=null) {
        withBody(new BodyDef(closure, source, section, values))
    }

    ProcessBuilder withBody(BodyDef body) {
        this.body = body
        return this
    }

    ProcessConfig getConfig() {
        return config
    }

    ProcessDef build() {
        if ( !body )
            throw new ScriptRuntimeException("Missing script in the specified process block -- make sure it terminates with the script string to be executed")
        return new ProcessDef(ownerScript, processName, body, config)
    }

}
