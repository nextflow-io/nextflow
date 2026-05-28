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

package nextflow.trace

import groovy.transform.CompileStatic

import java.util.concurrent.ConcurrentHashMap

/**
 * Generic key/value store for runtime workflow metadata that becomes
 * known after the workflow has started (e.g. identifiers assigned by
 * external schedulers, lineage URIs).
 *
 * Producers push values when they become authoritative; consumers
 * (telemetry observers) read an immutable snapshot. Neither side needs
 * to know the other; the bus is generic so new key/value pairs can be
 * added without changes to the transport layer.
 *
 * Storing {@code null} removes the key.
 */
@CompileStatic
class RuntimeMetadataStore {

    private final Map<String,Object> values = new ConcurrentHashMap<>()

    void put(String key, Object value) {
        if( !key )
            throw new IllegalArgumentException("Runtime metadata key cannot be null or empty")
        if( value == null )
            values.remove(key)
        else
            values.put(key, value)
    }

    Object get(String key) {
        return values.get(key)
    }

    Map<String,Object> snapshot() {
        return Collections.unmodifiableMap(new LinkedHashMap<String,Object>(values))
    }

    boolean isEmpty() {
        return values.isEmpty()
    }
}
