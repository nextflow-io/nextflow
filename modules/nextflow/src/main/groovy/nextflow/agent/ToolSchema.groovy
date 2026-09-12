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

/**
 * The pieces of tool-schema construction that every schema source shares: the JSON-schema
 * object envelope, the lenient type predicates, and the prose renderers that turn a declared
 * type into the words the model reads.
 *
 * <p>The sources themselves stay separate — {@link RecordSchema} reflects a record class,
 * {@link ProcessToolSchema} reads a typed process, {@link ModuleSpecToolSchema} reads a
 * sibling {@code meta.yml} and {@link ModuleMetadataToolSchema} reads the registry metadata —
 * because their input models and their type ladders genuinely differ. What they must NOT do is
 * spell the shared parts differently: every byte below is hashed into a tool's
 * {@code ToolDescriptor}, and from there into {@code AgentDef.toolsFingerprint}, the agent
 * {@code BodyDef.source} and the task hash. A cosmetic divergence here is a resume-cache
 * invalidation for every user wiring an affected module.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ToolSchema {

    /**
     * The JSON-schema object envelope, in a fresh mutable {@link LinkedHashMap}.
     *
     * <p>All four keys are written unconditionally and in this order — an empty
     * {@code properties} or {@code required} is emitted as an empty collection, never dropped.
     * A zero-input module and an all-optional record both produce empty accumulators, and the
     * {@code ls} tool declares a literal empty {@code required}; omitting the key would change
     * the schema JSON that {@code AgentDef.toolsFingerprint} hashes.
     *
     * <p>The returned map is fresh and mutable, so callers may hand it on as a nested fragment that
     * later code reads or extends. Only the OUTER map is fresh, though — {@code properties} and
     * {@code required} are stored by reference, not copied, so a caller that keeps hold of either
     * accumulator can still mutate the schema after the fact. Every caller today hands over an
     * accumulator it never touches again; a caller that wants to keep one must copy it itself.
     */
    static Map object(Map<String,Object> properties, List<String> required) {
        // a Groovy map literal IS a LinkedHashMap, so these four keys keep this order
        return [type: 'object', properties: properties, required: required, additionalProperties: false]
    }

    /**
     * The nf-core {@code meta.id} fragment: an open object carrying the single {@code id}
     * string property. There is no sub-schema for {@code meta} anywhere; {@code id} is the
     * convention mirrored from {@code CmdModuleView.inferNfCoreParam}.
     *
     * <p>WHEN it applies is the caller's business and the two callers disagree on purpose —
     * {@link ModuleSpecToolSchema} applies it to any {@code map} named {@code meta}, while
     * {@link ModuleMetadataToolSchema} gates it on the module being nf-core.
     */
    static Map metaIdFragment(String description) {
        return [
            type: 'object',
            description: description,
            properties: [id: [type: 'string', description: 'sample identifier']],
            additionalProperties: true ]
    }

    /** {@code file}/{@code path}/{@code directory} all denote an opaque path handle. */
    static boolean isFileType(String type) {
        return type == 'file' || type == 'path' || type == 'directory'
    }

    static boolean isIntegerType(String type) {
        return type == 'integer' || type == 'int' || type == 'long'
    }

    static boolean isNumberType(String type) {
        return type == 'float' || type == 'double' || type == 'number'
    }

    /**
     * Render a declared type as the phrase the model reads in a tool description.
     *
     * <p>{@code numberIsANumber} has NO default on purpose. The registry metadata maps
     * {@code float}/{@code double}/{@code number} onto their own rung; a {@code meta.yml}
     * spec does not, and lets them fall through to {@code 'a string'}. That divergence is
     * deliberate and out of scope here — it is on the description path, which is hashed into
     * the task key exactly as the schema path is, so unifying it silently re-runs every
     * affected agent on {@code -resume}. A default value is how the wrong caller would get the
     * wrong ladder without anyone choosing it.
     */
    static String describeKind(String type, boolean numberIsANumber) {
        final t = type?.toLowerCase()
        if( isFileType(t) )
            return 'a file path string'
        if( t == 'map' )
            return 'an object'
        if( isIntegerType(t) )
            return 'an integer'
        if( numberIsANumber && isNumberType(t) )
            return 'a number'
        if( t == 'boolean' )
            return 'a boolean'
        return 'a string'
    }

    /**
     * Render one output component as {@code `name` (kind)(description)}.
     *
     * <p>It takes already-extracted strings rather than a DTO so the two callers can keep their
     * own null handling: the registry items are null-navigated, the meta.yml components are not.
     */
    static String describeComponent(String name, String type, String description, boolean numberIsANumber) {
        final label = name ?: 'value'
        final kind = describeKind(type, numberIsANumber)
        final desc = description ? " (${description})" : ''
        return "`${label}` (${kind})${desc}".toString()
    }

}
