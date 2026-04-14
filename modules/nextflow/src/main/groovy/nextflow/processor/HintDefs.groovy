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
import groovy.util.logging.Slf4j

/**
 * Defines the global hint key registry and provides validation
 * for the {@code hints} process directive.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class HintDefs {

    /**
     * Global registry of known unprefixed hint keys and their expected value types.
     * Executor-prefixed keys (e.g. {@code seqera/machineRequirement.arch}) are
     * validated by the target executor, not this registry.
     */
    static final Map<String, Class> KNOWN_HINTS = [
        'consumableResource': String
    ]

    /**
     * Validates the given hints map against the global registry.
     * <p>
     * Unprefixed keys are checked against {@link #KNOWN_HINTS}:
     * <ul>
     *   <li>Unknown keys produce a warning with a "did you mean?" suggestion if a close match exists</li>
     *   <li>Value types are validated (must resolve to String or Integer)</li>
     * </ul>
     * Executor-prefixed keys (containing {@code /}) are skipped — they are validated by the target executor.
     *
     * @param hints the resolved hints map
     */
    static void validateHints(Map<String, Object> hints) {
        if( !hints )
            return

        for( Map.Entry<String, Object> entry : hints.entrySet() ) {
            final key = entry.key
            final value = entry.value

            // skip executor-prefixed keys — validated by the executor
            if( key.contains('/') )
                continue

            // validate value type
            if( value != null && !(value instanceof String) && !(value instanceof Integer) ) {
                throw new IllegalArgumentException("Invalid hint value type for key '${key}': expected String or Integer, got ${value.getClass().getName()}")
            }

            // validate key against registry
            if( !KNOWN_HINTS.containsKey(key) ) {
                final suggestions = KNOWN_HINTS.keySet().toList().closest(key)
                if( suggestions ) {
                    log.warn "Unknown process hint: '${key}' — did you mean '${suggestions.first()}'?"
                }
                else {
                    log.warn "Unknown process hint: '${key}'"
                }
            }
        }
    }

}
