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

package nextflow.processor

import groovy.transform.CompileStatic

/**
 * Validates the shape of the {@code hints} process directive.
 *
 * The core is intentionally agnostic about which hint keys are supported:
 * each executor validates the keys it recognizes (prefixed with its own
 * namespace, e.g. {@code awsbatch/...}, {@code seqera/...}). This class
 * only enforces that the map conforms to {@code Map<String,String>}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class HintDefs {

    /**
     * Validates the hint map structure. Does not check whether keys are
     * recognized — that is the responsibility of each executor.
     *
     * Rules:
     * <ul>
     *   <li>keys must be non-empty</li>
     *   <li>keys may contain at most one {@code /} separating the optional
     *       executor namespace from the hint name</li>
     *   <li>values must be {@code String} or {@code null}</li>
     * </ul>
     *
     * @param hints the hint map to validate (may be {@code null})
     * @throws IllegalArgumentException if the map is malformed
     */
    static void validateHints(Map<String, ?> hints) {
        if( !hints )
            return

        for( Map.Entry<String, ?> entry : hints.entrySet() ) {
            final key = entry.key
            final value = entry.value

            if( !key )
                throw new IllegalArgumentException("Process hint key cannot be null or empty")

            if( key.count('/') > 1 )
                throw new IllegalArgumentException("Invalid hint key '${key}': expected 'name' or 'executor/name'")

            if( value != null && !(value instanceof CharSequence) )
                throw new IllegalArgumentException("Invalid hint value for key '${key}': expected String, got ${value.getClass().getName()}")
        }
    }

}
