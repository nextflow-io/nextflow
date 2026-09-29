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
package nextflow.agent

import groovy.transform.CompileStatic
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpec.ModuleParam

/**
 * Derives a portable JSON-schema {@link Map} describing the inputs of a module
 * {@link ModuleSpec} (typically loaded from a sibling {@code meta.yml}), and a
 * human-readable description of its outputs.
 *
 * The schema is FLATTENED: a tuple input channel (e.g. {@code tuple(meta, reads)})
 * contributes one top-level property per component (so the LLM passes
 * {@code {"meta": {...}, "reads": "/abs/path"}} rather than a nested array), while
 * a scalar input channel contributes a single property. Each property carries the
 * {@code description} declared in the meta.yml so the LLM has the full context the
 * spec provides "for free".
 *
 * Type mapping is intentionally lenient (meta.yml types vary): {@code map} → object,
 * {@code file}/{@code path}/{@code directory} → string (path handle), numeric →
 * integer, {@code boolean} → boolean, everything else → string.
 *
 * langchain4j tools carry no output schema, so the LLM learns the output shape from
 * the tool description + the JSON it gets back; {@link #outputDescription} renders
 * that shape as prose to append to the tool description.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ModuleSpecToolSchema {

    /**
     * Build the flattened JSON-schema map for the module inputs. Every component of
     * every input channel becomes a top-level property keyed by its name; all are
     * marked {@code required} (basic meta.yml carries no optionality).
     *
     * @throws IllegalArgumentException if two flattened properties share a name
     */
    static Map inputSchema(ModuleSpec spec) {
        final properties = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()

        final inputs = spec.inputs ?: Collections.<ModuleParam>emptyList()
        for( final param : inputs ) {
            final components = param.isTuple() ? param.components : Collections.singletonList(param)
            for( final comp : components ) {
                final name = comp.name
                if( name == null || name.isEmpty() )
                    throw new IllegalArgumentException("Module spec `${spec.name}` declares an input with no name - cannot expose it as an agent tool property")
                if( properties.containsKey(name) )
                    throw new IllegalArgumentException("Module spec `${spec.name}` declares duplicate input name `${name}` across its input channels - agent tool properties must be unique")
                properties.put(name, fragmentFor(comp))
                required.add(name)
            }
        }

        return ToolSchema.object(properties, required)
    }

    /**
     * A human-readable description of the module outputs (names + component shapes)
     * to append to the tool description, so the LLM understands the JSON it gets back.
     * Files are described as absolute path strings (the opaque-path contract).
     */
    static String outputDescription(ModuleSpec spec) {
        final outputs = spec.outputs ?: Collections.<ModuleParam>emptyList()
        if( outputs.isEmpty() )
            return 'Returns a JSON object (no declared outputs).'

        final sb = new StringBuilder()
        sb.append('Returns a JSON object with the following output(s):')
        for( final param : outputs ) {
            final name = param.name ?: 'result'
            sb.append('\n- `').append(name).append('`: ')
            if( param.isTuple() ) {
                final parts = new ArrayList<String>()
                for( final comp : param.components )
                    parts.add(describeComponent(comp))
                sb.append('an object with ').append(parts.join(', '))
            }
            else {
                sb.append(describeComponent(param))
            }
        }
        sb.append('\nFile/path outputs are returned as absolute path strings (never file contents).')
        return sb.toString()
    }

    /**
     * A meta.yml component is never null-navigated: the spec loader materialises every
     * declared component. The {@code false} pins the number rung OFF - a {@code float} here
     * renders as {@code a string}, unlike the registry-sourced ladder.
     */
    private static String describeComponent(ModuleParam comp) {
        return ToolSchema.describeComponent(comp.name, comp.type, comp.description, false)
    }

    private static Map fragmentFor(ModuleParam param) {
        final type = param.type?.toLowerCase()
        final desc = param.description
        if( type == 'map' && 'meta'.equalsIgnoreCase(param.name) )
            return ToolSchema.metaIdFragment(desc ?: 'sample metadata')

        final Map<String,Object> fragment = new LinkedHashMap<String,Object>()
        if( type == 'map' ) {
            fragment.put('type', 'object')
            fragment.put('additionalProperties', true)
            if( desc )
                fragment.put('description', desc)
        }
        else if( ToolSchema.isFileType(type) ) {
            fragment.put('type', 'string')
            // make the path-handle contract explicit in the description
            fragment.put('description', desc ? "${desc} (file path)".toString() : 'file path')
        }
        else if( ToolSchema.isIntegerType(type) ) {
            fragment.put('type', 'integer')
            if( desc )
                fragment.put('description', desc)
        }
        else if( type == 'boolean' ) {
            fragment.put('type', 'boolean')
            if( desc )
                fragment.put('description', desc)
        }
        else {
            // val / string / unknown -> lenient string
            fragment.put('type', 'string')
            if( desc )
                fragment.put('description', desc)
        }
        return fragment
    }

}
