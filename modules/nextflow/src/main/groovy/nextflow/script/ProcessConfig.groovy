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

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskConfig
import nextflow.util.TestOnly

import static nextflow.util.CacheHelper.HashMode

/**
 * Holds the process configuration properties
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ProcessConfig implements Map<String,Object>, Cloneable {

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

    @Override
    Object getProperty( String name ) {

        switch( name ) {
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

    BaseScript getOwnerScript() { ownerScript }

    String getProcessName() { processName }

    TaskConfig createTaskConfig() {
        return new TaskConfig(configProperties)
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

    List<String> getSecret() {
        (List<String>) configProperties.get('secret') ?: Collections.<String>emptyList()
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

}
