/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.util

import java.lang.reflect.Field

import groovy.transform.CompileStatic
import io.seqera.config.MachineRequirementOpts
import nextflow.config.spec.ConfigOption

/**
 * Helper for processing {@code seqera/machineRequirement.*} hints from the
 * {@code hints} process directive and overlaying them onto
 * {@link MachineRequirementOpts} config-scope values.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class HintHelper {

    static final String PREFIX = 'seqera/'
    static final String MR_PREFIX = 'machineRequirement.'

    private static final List<Field> MR_FIELDS = Collections.unmodifiableList(
        MachineRequirementOpts.declaredFields
            .findAll { Field f -> f.isAnnotationPresent(ConfigOption) }
            .collect { Field f -> f.setAccessible(true); f } as List<Field>
    )

    static final Set<String> KNOWN_KEYS = Collections.unmodifiableSet(
        MR_FIELDS.collect { MR_PREFIX + it.name }.toSet()
    )

    private static final String SUPPORTED_KEYS_MSG =
        KNOWN_KEYS.collect { PREFIX + it }.sort().join(', ')

    /**
     * Extract hints consumed by the Seqera executor and validate them.
     *
     * <p>Both {@code seqera/}-prefixed and unprefixed keys that match one of the
     * {@link #KNOWN_KEYS} are returned with the prefix (if any) stripped. When the
     * same logical key appears both prefixed and unprefixed, the prefixed form
     * wins (executor-targeted hints override the general form).</p>
     *
     * <p>Foreign-namespaced keys (e.g. {@code awsbatch/...}) and unprefixed keys
     * that are not recognized are ignored — they may be targeted at another
     * executor. Unrecognized {@code seqera/}-prefixed keys raise an error, since
     * they were explicitly targeted at this executor.</p>
     *
     * @param hints the full hints map from task config
     * @return a map of known hint names (no prefix) to values
     */
    static Map<String, Object> extractSeqeraHints(Map<String, Object> hints) {
        if( !hints )
            return Collections.emptyMap()

        final unprefixed = new LinkedHashMap<String, Object>()
        final prefixed = new LinkedHashMap<String, Object>()
        for( Map.Entry<String, Object> entry : hints.entrySet() ) {
            final key = entry.key
            if( !key )
                continue

            if( key.startsWith(PREFIX) ) {
                final stripped = key.substring(PREFIX.length())
                if( !KNOWN_KEYS.contains(stripped) )
                    throw new IllegalArgumentException("Unknown Seqera Platform hint: '${key}' — supported keys are: ${SUPPORTED_KEYS_MSG}")
                prefixed.put(stripped, entry.value)
            }
            else if( !key.contains('/') && KNOWN_KEYS.contains(key) ) {
                unprefixed.put(key, entry.value)
            }
        }

        unprefixed.putAll(prefixed)
        return unprefixed
    }

    /**
     * Overlay {@code machineRequirement.*} hints onto existing config-scope
     * {@link MachineRequirementOpts}. Hint values take precedence over
     * config-scope values.
     *
     * @param baseOpts the config-scope machine requirement options
     * @param hints the full hints map from task config
     * @return a new {@link MachineRequirementOpts} with hints overlaid
     */
    static MachineRequirementOpts overlayHints(MachineRequirementOpts baseOpts, Map<String, Object> hints) {
        final seqeraHints = extractSeqeraHints(hints)
        if( !seqeraHints )
            return baseOpts

        final Map<String, Object> merged = new LinkedHashMap<>()
        for( final field : MR_FIELDS ) {
            final value = field.get(baseOpts)
            if( value != null )
                merged.put(field.name, value)
        }

        for( Map.Entry<String, Object> entry : seqeraHints.entrySet() ) {
            final fieldName = entry.key.substring(MR_PREFIX.length())
            final value = entry.value
            if( value == null ) {
                merged.remove(fieldName)
                continue
            }
            merged.put(fieldName, value)
        }

        return new MachineRequirementOpts(merged)
    }

}
