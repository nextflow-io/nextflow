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

import groovy.transform.CompileStatic
import io.seqera.config.MachineRequirementOpts

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

    /**
     * Known {@code seqera/machineRequirement.*} hint keys (after stripping the {@code seqera/} prefix).
     */
    static final Set<String> KNOWN_KEYS = Set.of(
        'machineRequirement.arch',
        'machineRequirement.provisioning',
        'machineRequirement.maxSpotAttempts',
        'machineRequirement.machineTypes',
        'machineRequirement.diskType',
        'machineRequirement.diskThroughputMiBps',
        'machineRequirement.diskIops',
        'machineRequirement.diskEncrypted',
        'machineRequirement.diskAllocation',
        'machineRequirement.diskMountPath',
        'machineRequirement.diskSize',
        'machineRequirement.capacityMode'
    )

    /**
     * Extract {@code seqera/}-prefixed hints from the hints map and validate them.
     * Raises an error for unrecognized {@code seqera/} keys.
     *
     * @param hints the full hints map from task config
     * @return a map of hint keys (with {@code seqera/} prefix stripped) to values
     */
    static Map<String, String> extractSeqeraHints(Map<String, String> hints) {
        if( !hints )
            return Collections.emptyMap()

        final result = new LinkedHashMap<String, String>()
        for( Map.Entry<String, String> entry : hints.entrySet() ) {
            final key = entry.key
            if( !key.startsWith(PREFIX) )
                continue

            final stripped = key.substring(PREFIX.length())
            if( !KNOWN_KEYS.contains(stripped) ) {
                throw new IllegalArgumentException("Unknown Seqera Platform hint: '${key}' — supported keys are: ${KNOWN_KEYS.collect { PREFIX + it }.sort().join(', ')}")
            }
            result.put(stripped, entry.value)
        }
        return result
    }

    /**
     * Overlay {@code seqera/machineRequirement.*} hints onto existing config-scope
     * {@link MachineRequirementOpts}. Hint values take precedence over config-scope values.
     *
     * @param baseOpts the config-scope machine requirement options
     * @param hints the full hints map from task config
     * @return a new {@link MachineRequirementOpts} with hints overlaid
     */
    static MachineRequirementOpts overlayHints(MachineRequirementOpts baseOpts, Map<String, String> hints) {
        final seqeraHints = extractSeqeraHints(hints)
        if( !seqeraHints )
            return baseOpts

        // build a map from current opts, then overlay hints
        final Map<String, Object> merged = new LinkedHashMap<>()

        // copy existing config-scope values
        if( baseOpts.arch ) merged.arch = baseOpts.arch
        if( baseOpts.provisioning ) merged.provisioning = baseOpts.provisioning
        if( baseOpts.maxSpotAttempts != null ) merged.maxSpotAttempts = baseOpts.maxSpotAttempts
        if( baseOpts.machineTypes ) merged.machineTypes = baseOpts.machineTypes
        if( baseOpts.diskType ) merged.diskType = baseOpts.diskType
        if( baseOpts.diskThroughputMiBps != null ) merged.diskThroughputMiBps = baseOpts.diskThroughputMiBps
        if( baseOpts.diskIops != null ) merged.diskIops = baseOpts.diskIops
        if( baseOpts.diskEncrypted != null ) merged.diskEncrypted = baseOpts.diskEncrypted
        if( baseOpts.diskAllocation ) merged.diskAllocation = baseOpts.diskAllocation
        if( baseOpts.diskMountPath ) merged.diskMountPath = baseOpts.diskMountPath
        if( baseOpts.diskSize ) merged.diskSize = baseOpts.diskSize.toString()
        if( baseOpts.capacityMode ) merged.capacityMode = baseOpts.capacityMode

        // overlay hints — strip the machineRequirement. prefix to get field names
        for( Map.Entry<String, String> entry : seqeraHints.entrySet() ) {
            final fieldName = entry.key.substring(MR_PREFIX.length())
            final value = entry.value
            // handle special type conversions
            switch( fieldName ) {
                case 'machineTypes':
                    // comma-separated string → List<String>
                    merged.machineTypes = value.split(',').collect { it.trim() }
                    break
                case 'diskEncrypted':
                    merged.diskEncrypted = Boolean.parseBoolean(value)
                    break
                default:
                    merged.put(fieldName, value)
                    break
            }
        }

        return new MachineRequirementOpts(merged)
    }

}
