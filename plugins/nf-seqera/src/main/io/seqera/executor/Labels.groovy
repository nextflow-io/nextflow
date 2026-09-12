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

package io.seqera.executor

import groovy.transform.CompileStatic

/**
 * Helper class to manage run labels.
 *
 * Builds the labels map from the resource labels attached to the run -- the auto-derived
 * workflow metadata labels, see {@link nextflow.platform.AutoLabels}, and the config-level
 * {@code process.resourceLabels}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Labels {

    /**
     * Merge the two label maps, coercing keys and values to string via {@link #toStringMap}.
     * These are the labels auto-derived from the workflow metadata ({@code base}), overlaid by
     * the config-level {@code process.resourceLabels} ({@code overlay}) which win on a key collision.
     */
    static Map<String,String> merge(Map<String,?> base, Map<String,?> overlay) {
        final result = new LinkedHashMap<String,String>(toStringMap(base))
        result.putAll(toStringMap(overlay))
        return result
    }

    /**
     * Coerce arbitrary map values to strings via {@link String#valueOf}.
     * Returns an empty map for null/empty input. Throws
     * {@link IllegalArgumentException} when the value is not a {@link Map},
     * to surface a clear error when {@code process.resourceLabels} is
     * misconfigured (e.g. as a list).
     */
    static Map<String,String> toStringMap(Object value) {
        if( value == null )
            return Collections.<String,String>emptyMap()
        if( value !instanceof Map )
            throw new IllegalArgumentException("Invalid value for 'resourceLabels' directive - expected a map of key/value pairs, got '${value.getClass().getName()}'")
        final map = (Map<?,?>) value
        if( map.isEmpty() )
            return Collections.<String,String>emptyMap()
        final result = new LinkedHashMap<String,String>(map.size())
        for( Map.Entry<?,?> entry : map.entrySet() )
            result.put(entry.key.toString(), String.valueOf(entry.value))
        return result
    }

    /**
     * Return the entries of {@code task} that are missing from {@code run}
     * or have a different value. Returns {@code null} if the resulting
     * map would be empty (so callers can omit the field).
     */
    static Map<String,String> delta(Map<String,String> task, Map<String,String> run) {
        if( !task ) return null
        final result = new LinkedHashMap<String,String>()
        for( Map.Entry<String,String> entry : task.entrySet() ) {
            final k = entry.key
            final v = entry.value
            if( run == null || !run.containsKey(k) || run.get(k) != v )
                result.put(k, v)
        }
        return result.isEmpty() ? null : result
    }
}
