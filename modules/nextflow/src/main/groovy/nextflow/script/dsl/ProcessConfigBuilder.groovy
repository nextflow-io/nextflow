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
import nextflow.exception.ConfigParseException
import nextflow.script.ProcessConfig

/**
 * Builder for {@link ProcessConfig}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@TypeChecked
class ProcessConfigBuilder extends ProcessBuilder {

    ProcessConfigBuilder(ProcessConfig config) {
        super(config)
    }

    /**
     * Apply process config settings from the config file to a process.
     *
     * @param configProcessScope
     * @param baseName
     * @param simpleName
     * @param fullyQualifiedName
     */
    void applyConfig(Map configProcessScope, String baseName, String simpleName, String fullyQualifiedName) {
        // -- apply settings defined in the config object using the`withLabel:` syntax
        final processLabels = config.getLabels() ?: ['']
        applyConfigSelectorWithLabels(configProcessScope, processLabels)

        // -- apply settings defined in the config file using the process base name
        applyConfigSelectorWithName(configProcessScope, baseName)

        // -- apply settings defined in the config file using the process simple name
        if( simpleName && simpleName!=baseName )
            applyConfigSelectorWithName(configProcessScope, simpleName)

        // -- apply settings defined in the config file using the process fully qualified name (ie. with the execution scope)
        if( fullyQualifiedName && (fullyQualifiedName!=simpleName || fullyQualifiedName!=baseName) )
            applyConfigSelectorWithName(configProcessScope, fullyQualifiedName)

        // -- apply defaults
        applyConfigDefaults(configProcessScope)

        // -- check for conflicting settings
        if( config.scratch && config.stageInMode == 'rellink' ) {
            log.warn("Directives `scratch` and `stageInMode=rellink` conflict with each other -- Enforcing default stageInMode for process `$simpleName`")
            config.remove('stageInMode')
        }
    }

    /**
     * Apply the config settings in a label selector, for example:
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
     * @param labels
     */
    protected void applyConfigSelectorWithLabels(Map<String,?> configDirectives, List<String> labels) {
        final prefix = 'withLabel:'
        for( String rule : configDirectives.keySet() ) {
            if( !rule.startsWith(prefix) )
                continue
            final pattern = rule.substring(prefix.size()).trim()
            if( !matchesLabels(labels, pattern) )
                continue

            log.debug "Config settings `$rule` matches labels `${labels.join(',')}` for process with name $processName"
            final settings = configDirectives.get(rule)
            if( settings instanceof Map ) {
                applyConfigSettings(settings)
            }
            else if( settings != null ) {
                throw new ConfigParseException("Unknown config settings for process labeled ${labels.join(',')} -- settings=$settings ")
            }
        }
    }

    /**
     * Determine whether any of the given labels match the
     * given process withLabel selector.
     *
     * @param labels
     * @param pattern
     */
    static boolean matchesLabels(List<String> labels, String pattern) {
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

    /**
     * Apply the config settings in a name selector, for example:
     *
     * ```
     * process {
     *     withName: foo {
     *         cpus = 1
     *         memory = 2.gb
     *     }
     * }
     * ```
     *
     * @param configDirectives
     * @param target
     */
    protected void applyConfigSelectorWithName(Map<String,?> configDirectives, String target) {
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

    /**
     * Determine whether the given process name matches the
     * given process withName selector.
     *
     * @param name
     * @param pattern
     */
    static boolean matchesSelector(String name, String pattern) {
        final isNegated = pattern.startsWith('!')
        if( isNegated )
            pattern = pattern.substring(1).trim()
        return Pattern.compile(pattern).matcher(name).matches() ^ isNegated
    }

    /**
     * Apply config settings to a process.
     *
     * @param settings
     */
    protected void applyConfigSettings(Map<String,?> settings) {
        if( !settings )
            return

        for( final entry : settings ) {
            if( entry.key.startsWith("withLabel:") || entry.key.startsWith("withName:"))
                continue

            if( !DIRECTIVES.contains(entry.key) )
                log.warn "Unknown directive `$entry.key` for process `$processName`"

            if( entry.key == 'params' ) // <-- patch issue #242
                continue

            if( entry.key == 'ext' ) {
                final ext = config.getProperty('ext')
                if( ext instanceof Map ) {
                    // update missing 'ext' properties found in 'process' scope
                    ext.putAll(entry.value as Map)
                }
                continue
            }

            putWithRepeat(entry.key, entry.value)
        }
    }

    /**
     * Apply the global settings in the process config scope to a process.
     *
     * @param defaults
     *      A map object representing the setting to be applied to the process
     *      (provided it does not already define a different value for
     *      the same config setting).
     */
    protected void applyConfigDefaults( Map defaults ) {
        for( String key : defaults.keySet() ) {
            if( key == 'params' )
                continue
            final value = defaults.get(key)
            final current = config.getProperty(key)
            if( key == 'ext' ) {
                if( value instanceof Map && current instanceof Map ) {
                    final ext = current as Map
                    value.each { k,v -> if(!ext.containsKey(k)) ext.put(k,v) }
                }
            }
            else if( !config.containsKey(key) || (ProcessConfig.DEFAULT_CONFIG.containsKey(key) && current==ProcessConfig.DEFAULT_CONFIG.get(key)) ) {
                putWithRepeat(key, value)
            }
        }
    }

    private static final List<String> REPEATABLE_DIRECTIVES = ['label','module','pod','publishDir']

    protected void putWithRepeat( String name, Object value ) {
        if( name in REPEATABLE_DIRECTIVES ) {
            config.remove(name)
            this.metaClass.invokeMethod(this, name, value)
        }
        else {
            config.put(name, value)
        }
    }

}
