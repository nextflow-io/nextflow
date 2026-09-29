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
import io.seqera.npr.api.schema.v1.ModuleChannel
import io.seqera.npr.api.schema.v1.ModuleChannelItem
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.api.schema.v1.ModuleTool

/**
 * Derives a portable JSON-schema {@link Map} describing the inputs of a module from the
 * registry {@link ModuleMetadata} (fetched from {@code GET /api/v1/modules/{name}}), and a
 * human-readable tool description from the module + tool + output metadata.
 *
 * This is the registry-sourced counterpart of {@link ModuleSpecToolSchema} (which works from
 * a sibling {@code meta.yml} {@link nextflow.module.ModuleSpec}). It exists because the registry
 * metadata is richer (per-field descriptions, patterns, enums, tool homepages/documentation) and
 * is the canonical, always-available source when the module is resolved from the registry.
 *
 * The schema is FLATTENED: each {@link ModuleChannel} item contributes one top-level property
 * keyed by its name (so the LLM passes {@code {"meta": {...}, "reads": "/abs/path"}} rather than
 * a nested array), carrying its {@code description} (plus any {@code pattern}/{@code enum}).
 *
 * <p><b>nf-core {@code meta.id} convention:</b> for an nf-core module, a {@code map} item named
 * {@code meta} is exposed as a nested object schema with an {@code id} string property — mirroring
 * the hardcoded {@code --meta.id <ID>} convention in {@code CmdModuleView.inferNfCoreParam}. There
 * is no sub-schema for {@code meta} in the registry; the {@code id} property is the convention.
 *
 * <p>All getters on the npr DTOs can return {@code null}; everything is null-guarded.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ModuleMetadataToolSchema {

    /**
     * Build the flattened JSON-schema map for the module inputs from the registry metadata.
     * Every item of every input channel becomes a top-level property keyed by its name; all
     * are marked {@code required}.
     *
     * @param metadata the registry module metadata
     * @param nfCore   whether the module reference scope is {@code nf-core} (enables the
     *                 {@code meta.id} convention)
     * @throws IllegalArgumentException if an input item has no name or two flattened
     *         properties share a name
     */
    static Map inputSchema(ModuleMetadata metadata, boolean nfCore) {
        final properties = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()

        final inputs = metadata?.getInput() ?: Collections.<ModuleChannel>emptyList()
        for( final ModuleChannel channel : inputs ) {
            if( channel == null )
                continue
            final items = channel.getItems() ?: Collections.<ModuleChannelItem>emptyList()
            for( final ModuleChannelItem item : items ) {
                if( item == null )
                    continue
                final name = item.getName()
                if( name == null || name.isEmpty() )
                    throw new IllegalArgumentException("Module metadata declares an input item with no name - cannot expose it as an agent tool property")
                if( properties.containsKey(name) )
                    throw new IllegalArgumentException("Module metadata declares duplicate input name `${name}` across its input channels - agent tool properties must be unique")
                properties.put(name, fragmentFor(item, nfCore))
                required.add(name)
            }
        }

        return ToolSchema.object(properties, required)
    }

    /**
     * The flattened input property names derived from the metadata, in declaration order.
     * Used by the bridge to cross-check against the executable {@code ModuleSpec} (drift guard).
     */
    static List<String> inputPropertyNames(ModuleMetadata metadata) {
        final names = new ArrayList<String>()
        final inputs = metadata?.getInput() ?: Collections.<ModuleChannel>emptyList()
        for( final ModuleChannel channel : inputs ) {
            if( channel == null )
                continue
            final items = channel.getItems() ?: Collections.<ModuleChannelItem>emptyList()
            for( final ModuleChannelItem item : items ) {
                final name = item?.getName()
                if( name )
                    names.add(name)
            }
        }
        return names
    }

    /**
     * Build a JSON-schema fragment for a single input item. The mapping mirrors
     * {@link ModuleSpecToolSchema}, folding in {@code description}, {@code pattern} (appended to the
     * description) and {@code enum}. The nf-core {@code meta} map item gets a nested {@code id}
     * property (the {@code meta.id} convention).
     */
    private static Map fragmentFor(ModuleChannelItem item, boolean nfCore) {
        return buildFragment(
            item.getName(),
            item.getType(),
            item.getDescription(),
            item.getPattern(),
            item.get_enum(),
            nfCore )
    }

    /**
     * The per-item mapping seam (primitive inputs) - keeps the npr-DTO coupling at the edge and
     * lets the mapping logic be unit-tested without constructing a {@link ModuleMetadata}.
     */
    protected static Map buildFragment(String name, String rawType, String description, String pattern, List<String> enumValues, boolean nfCore) {
        final type = rawType?.toLowerCase()

        // -- nf-core `meta.id` convention: a `map` item named `meta` for an nf-core module
        //    is exposed as a nested object with an `id` string property (no sub-schema exists
        //    in the registry; `id` is the hardcoded convention mirrored from CmdModuleView).
        if( nfCore && type == 'map' && 'meta'.equalsIgnoreCase(name) )
            return ToolSchema.metaIdFragment(composeDescription(description ?: 'sample metadata', pattern))

        final Map<String,Object> fragment = new LinkedHashMap<String,Object>()
        if( type == 'map' ) {
            fragment.put('type', 'object')
            final desc = composeDescription(description, pattern)
            if( desc )
                fragment.put('description', desc)
            fragment.put('additionalProperties', true)
            return fragment
        }

        if( ToolSchema.isFileType(type) ) {
            fragment.put('type', 'string')
            final base = description ? "${description} (file path)".toString() : 'file path'
            fragment.put('description', composeDescription(base, pattern))
        }
        else if( ToolSchema.isIntegerType(type) ) {
            fragment.put('type', 'integer')
            putDescriptionIfPresent(fragment, description, pattern)
        }
        else if( ToolSchema.isNumberType(type) ) {
            fragment.put('type', 'number')
            putDescriptionIfPresent(fragment, description, pattern)
        }
        else if( type == 'boolean' ) {
            fragment.put('type', 'boolean')
            putDescriptionIfPresent(fragment, description, pattern)
        }
        else {
            // string / val / unknown -> lenient string
            fragment.put('type', 'string')
            putDescriptionIfPresent(fragment, description, pattern)
        }

        if( enumValues != null && !enumValues.isEmpty() )
            fragment.put('enum', new ArrayList<String>(enumValues))

        return fragment
    }

    private static void putDescriptionIfPresent(Map<String,Object> fragment, String description, String pattern) {
        final desc = composeDescription(description, pattern)
        if( desc )
            fragment.put('description', desc)
    }

    /** Append a {@code (pattern: ...)} hint to the description when a pattern is present. */
    private static String composeDescription(String description, String pattern) {
        if( !pattern )
            return description
        return description
            ? "${description} (pattern: ${pattern})".toString()
            : "pattern: ${pattern}".toString()
    }

    /**
     * A human-readable tool description: the module description, the wrapped tools (name +
     * homepage/documentation URIs) and the output shape (the {@code output} Map keys are the
     * emit names, each with its component items). Everything is null-guarded.
     */
    static String description(ModuleMetadata metadata) {
        final sb = new StringBuilder()
        final desc = metadata?.getDescription()
        sb.append(desc ?: 'module tool')

        final tools = metadata?.getTools()
        if( tools ) {
            for( final ModuleTool tool : tools ) {
                if( tool == null )
                    continue
                final name = tool.getName()
                if( !name )
                    continue
                sb.append('\nTool `').append(name)
                if( tool.getVersion() )
                    sb.append("` v").append(tool.getVersion())
                else
                    sb.append('`')
                final links = new ArrayList<String>()
                if( tool.getHomepage() )
                    links.add("homepage: ${tool.getHomepage()}".toString())
                if( tool.getDocumentation() )
                    links.add("documentation: ${tool.getDocumentation()}".toString())
                if( links )
                    sb.append(' (').append(links.join(', ')).append(')')
            }
        }

        sb.append('\n\n').append(outputDescription(metadata))
        return sb.toString()
    }

    /**
     * A human-readable description of the module outputs (emit names + component shapes) to
     * append to the tool description. The registry {@code output} is a Map keyed by emit name.
     */
    static String outputDescription(ModuleMetadata metadata) {
        final outputs = metadata?.getOutput() ?: Collections.<String,ModuleChannel>emptyMap()
        if( outputs.isEmpty() )
            return 'Returns a JSON object (no declared outputs).'

        final sb = new StringBuilder()
        sb.append('Returns a JSON object with the following output(s):')
        for( final Map.Entry<String,ModuleChannel> entry : outputs.entrySet() ) {
            final emit = entry.key ?: 'result'
            final channel = entry.value
            sb.append('\n- `').append(emit).append('`: ')
            final items = channel?.getItems() ?: Collections.<ModuleChannelItem>emptyList()
            if( items.isEmpty() ) {
                sb.append('a value')
                continue
            }
            if( items.size() > 1 || Boolean.TRUE.equals(channel?.getTuple()) ) {
                final parts = new ArrayList<String>()
                for( final ModuleChannelItem item : items )
                    parts.add(describeItem(item))
                sb.append('an object with ').append(parts.join(', '))
            }
            else {
                sb.append(describeItem(items[0]))
            }
        }
        sb.append('\nFile/path outputs are returned as absolute path strings (never file contents).')
        return sb.toString()
    }

    /**
     * Every registry DTO getter can return {@code null}, so the item is null-navigated HERE and
     * the shared renderer only ever sees extracted strings. The {@code true} pins the number
     * rung ON - the registry declares {@code float}/{@code double}/{@code number} types that a
     * meta.yml does not.
     */
    private static String describeItem(ModuleChannelItem item) {
        return ToolSchema.describeComponent(item?.getName(), item?.getType(), item?.getDescription(), true)
    }

}
