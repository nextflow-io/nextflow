/*
 * Copyright 2013-2025, Seqera Labs
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

import groovy.transform.TypeChecked
import groovy.util.logging.Slf4j
import nextflow.exception.IllegalConfigException
import nextflow.exception.IllegalDirectiveException
import nextflow.exception.ScriptRuntimeException
import nextflow.processor.ConfigList
import nextflow.processor.ErrorStrategy
import nextflow.script.BaseScript
import nextflow.script.BodyDef
import nextflow.script.ProcessConfig
import nextflow.script.ProcessDef

/**
 * Builder for {@link ProcessDef}.
 *
 * @see nextflow.script.dsl.ProcessDsl.DirectiveDsl
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@TypeChecked
class ProcessBuilder {

    static final List<String> DIRECTIVES = [
            'accelerator',
            'afterScript',
            'arch',
            'array',
            'beforeScript',
            'cache',
            'clusterOptions',
            'conda',
            'consumableResources',
            'container',
            'containerOptions',
            'cpus',
            'debug',
            'disk',
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
            'resourceLimits',
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

    ProcessBuilder(ProcessConfig config) {
        this.ownerScript = config.getOwnerScript()
        this.processName = config.getProcessName()
        this.config = config
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
        if( name == 'when' || name == 'stub' )
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

    void accelerator( Map params, value ) {
        if( value instanceof Number ) {
            if( params.limit==null )
                params.limit=value
            else if( params.request==null )
                params.request=value
        }
        else if( value != null ) {
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
        }
        accelerator(params)
    }

    void accelerator( value ) {
        if( value instanceof Number )
            config.put('accelerator', [limit: value])
        else if( value instanceof Map || value instanceof Closure )
            config.put('accelerator', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `accelerator` directive value: $value [${value.getClass().getName()}]")
    }

    void arch( Map params, value ) {
        if( value instanceof String ) {
            if( params.name==null )
                params.name=value
        }
        else if( value != null ) {
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
        }
        arch(params)
    }

    void arch( value ) {
        if( value instanceof String )
            config.put('arch', [name: value])
        else if( value instanceof Map || value instanceof Closure )
            config.put('arch', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `arch` directive value: $value [${value.getClass().getName()}]")
    }

    void consumableResources( value ) {
        if( value instanceof List || value instanceof Closure )
            config.put('consumableResources', value)
        else if( value != null )
            throw new IllegalArgumentException("Not a valid `consumableResources` directive value: $value [${value.getClass().getName()}] - expected a list of maps with 'type' and 'value' keys")
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
    void disk( Map opts, value ) {
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
    void errorStrategy( strategy ) {
        if( strategy instanceof CharSequence && !ErrorStrategy.isValid(strategy) )
            throw new IllegalArgumentException("Unknown error strategy '${strategy}' ― Available strategies are: ${ErrorStrategy.values().join(',').toLowerCase()}")

        config.put('errorStrategy', strategy)
    }

    /**
     * Implements the {@code label} directive.
     *
     * This directive can be specified (invoked) more than once in
     * the process definition.
     *
     * @param value
     */
    void label(String value) {
        if( !value ) return

        // -- check that label has a valid syntax
        if( !isValidLabel(value) )
            throw new IllegalConfigException("Not a valid process label: $value -- Label must consist of alphanumeric characters or '_', must start with an alphabetic character and must end with an alphanumeric character")

        // -- get the current label, it must be a list
        def allLabels = (List)config.get('label')
        if( !allLabels ) {
            allLabels = new ConfigList()
            config.put('label', allLabels)
        }

        // -- avoid duplicates
        if( !allLabels.contains(value) )
            allLabels.add(value)
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
    void module( value ) {
        if( !value ) return

        def result = (List)config.module
        if( result == null ) {
            result = new ConfigList()
            config.put('module', result)
        }

        if( value instanceof List )
            result.addAll(value)
        else
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
            allOptions = new ConfigList()
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
            dirs = new ConfigList()
            config.put('publishDir', dirs)
        }

        dirs.add(params)
    }

    /**
     * Implements the {@code publishDir} directive as a path with named parameters, eg:
     *
     *     publishDir '/some/dir', mode: 'copy'
     *
     * @param params map of publish options
     * @param path   String | Closure<String>
     */
    void publishDir(Map params, path) {
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

    private static final List<String> VALID_RESOURCE_LIMITS = List.of('cpus', 'memory', 'disk', 'time')

    /**
     * Implements the {@code resourceLimits} directive.
     *
     * @param entries
     */
    void resourceLimits(Map entries) {
        for( entry in entries ) {
            if( entry.key !in VALID_RESOURCE_LIMITS )
                throw new IllegalArgumentException("Not a valid directive in `resourceLimits`: $entry.key")
        }

        config.put('resourceLimits', entries)
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
            allSecrets = new ConfigList()
            config.put('secret', allSecrets)
        }

        // -- avoid duplicates
        if( !allSecrets.contains(name) )
            allSecrets.add(name)
    }

    /// SCRIPT

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
        return new ProcessDef(ownerScript, processName, config, body)
    }

}
