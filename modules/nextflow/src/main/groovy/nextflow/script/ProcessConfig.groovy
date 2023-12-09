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

package nextflow.script

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskConfig
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
            maxRetries: 0,
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
     * List of process input definitions
     */
    private ProcessInputs inputs

    /**
     * List of process output definitions
     */
    private ProcessOutputs outputs

    protected ProcessConfig( BaseScript script ) {
        ownerScript = script
        configProperties = new LinkedHashMap()
        configProperties.putAll( DEFAULT_CONFIG )
    }

    ProcessConfig( BaseScript script, String name ) {
        this(script)
        this.processName = name
    }

    /* Only for testing purpose */
    @PackageScope
    ProcessConfig( Map delegate ) {
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

    ProcessConfig setInputs(ProcessInputs inputs) {
        this.inputs = inputs
        return this
    }

    ProcessConfig setOutputs(ProcessOutputs outputs) {
        this.outputs = outputs
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
                if( !configProperties.containsKey(name) ) {
                    configProperties.put(name, new HashMap())
                }
                return configProperties.get(name)

            default:
                if( configProperties.containsKey(name) )
                    return configProperties.get(name)
                else
                    return null
        }

    }

    @PackageScope
    BaseScript getOwnerScript() { ownerScript }

    @PackageScope
    String getProcessName() { processName }

    TaskConfig createTaskConfig() {
        return new TaskConfig(configProperties)
    }

    ProcessInputs getInputs() {
        inputs
    }

    ProcessOutputs getOutputs() {
        outputs
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

}
