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

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import nextflow.exception.ScriptRuntimeException

/**
 * A single entry of the agent {@code tools} directive, parsed and validated against the
 * declaration grammar: {@code family[:segment]*:name}.
 *
 * <p>This type is purely <b>syntactic</b>. It knows nothing about which families exist or
 * which tools they contain — that is {@link ToolRefResolver}'s job — so a ref that parses
 * here may still be rejected later as unknown or as matching nothing. Keeping the two apart
 * is what lets the grammar be unit-tested without a session, and what makes the five
 * zero-match failures distinguishable in the error message the user sees.
 *
 * <p>The rules enforced here, in the order they are checked:
 * <ul>
 *   <li><b>G6</b> — a leading {@code !} (negation) is rejected outright: there is no exclude
 *       operator, the safety boundary sits on the {@code shell:} family line instead;</li>
 *   <li><b>G5</b> — a glob must be anchored to a family, so a bare {@code *} is rejected;</li>
 *   <li><b>G1/G2</b> — at least two colon-separated segments, none of them empty. The value is
 *       taken <b>verbatim</b>: it is never trimmed, so {@code 'fs:read '} is an error rather
 *       than a silently repaired {@code fs:read};</li>
 *   <li><b>G4</b> — {@code *} may appear only in the last segment; the family and every
 *       intermediate segment must name a node exactly.</li>
 * </ul>
 *
 * <p>Matching is <b>case-sensitive</b> in every segment (see {@link #matches}), so
 * {@code nf:module_run:samtools_*} does not select a process named {@code SAMTOOLS_SORT}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'ref', includeFields = true)
class ToolRef {

    /** The one wildcard character the grammar defines; legal in the last segment only (G4). */
    static final String GLOB = '*'

    /**
     * Legal characters of a non-terminal segment. Deliberately narrower than a Nextflow
     * identifier: the terminal segment becomes a wire name, which OpenAI restricts to
     * {@code [a-zA-Z0-9_-]} (§4), and there is no reason for the declaration namespace to
     * admit characters the wire namespace cannot carry.
     */
    private static final Pattern PLAIN_SEGMENT = ~/[A-Za-z0-9_-]+/

    /** Same, plus the wildcard — the terminal segment only. */
    private static final Pattern GLOB_SEGMENT = ~/[A-Za-z0-9_*-]+/

    private final String ref

    private final List<String> segments

    private ToolRef(String ref, List<String> segments) {
        this.ref = ref
        this.segments = segments
    }

    /** The entry exactly as declared. */
    String getRef() { ref }

    /** The colon-separated segments; always at least two, none empty. */
    List<String> getSegments() { segments }

    /** The family this ref is anchored to, i.e. the first segment. */
    String getFamily() { segments[0] }

    /** The terminal segment — the only one that may carry a glob. */
    private String getLeafPattern() { segments[segments.size() - 1] }

    /** Whether the terminal segment is a pattern rather than an exact name. */
    boolean isGlobbed() { leafPattern.contains(GLOB) }

    @Override
    String toString() { ref }

    /**
     * Parse and validate a declared entry. The argument is the <b>resolved</b> value of the
     * directive entry ({@code entry?.toString()}), not the source literal, so an interpolated
     * entry is validated per invocation (G1).
     *
     * @throws ScriptRuntimeException when the value is not a well-formed ref
     */
    static ToolRef parse(String value) {
        if( !value )
            throw invalid(value, 'a tool reference cannot be empty')
        // -- G6: there is no exclude operator. Caught before the shape checks so the user is
        //    told the operator does not exist rather than that `!fs` is an odd family name.
        if( value.startsWith('!') )
            throw invalid(value, 'exclusions are not supported - drop the leading `!` and declare only the tools the agent may use')

        final String[] parts = value.split(':', -1)
        if( parts.length < 2 ) {
            // an unanchored glob is its own failure (G5): with no family it means "every tool
            // that exists", including ones a later release adds
            if( value.contains(GLOB) )
                throw invalid(value, 'a glob must be anchored to a tool family, e.g. `fs:*` or `nf:module_run:*`')
            throw invalid(value, 'a tool reference must be namespaced as `family[:group]:name`, e.g. `nf:module_run:MY_PROCESS`, `fs:*` or `shell:bash` - to expose a module, `include` it and name its process')
        }

        for( int i = 0; i < parts.length; i++ ) {
            final seg = parts[i]
            final last = i == parts.length - 1
            if( !seg )
                throw invalid(value, "segment ${i + 1} is empty - segments are separated by a single colon" as String)
            // -- G4: the family and every intermediate segment must name a node exactly, so the
            //    ref always states which family (and which group within it) it selects from
            if( !last && seg.contains(GLOB) )
                throw invalid(value, "`${seg}` cannot contain a glob - only the last segment may" as String)
            final ok = last
                ? GLOB_SEGMENT.matcher(seg).matches()
                : PLAIN_SEGMENT.matcher(seg).matches()
            if( !ok )
                throw invalid(value, "`${seg}` is not a legal segment - only letters, digits, `_` and `-` are allowed (plus `*` in the last segment), and an entry is never trimmed" as String)
        }

        return new ToolRef(value, List.<String>copyOf(Arrays.asList(parts)))
    }

    private static ScriptRuntimeException invalid(String value, String reason) {
        return new ScriptRuntimeException("Invalid tool reference `${value}` - ${reason}")
    }

    /**
     * Match one segment pattern against one node name, case-sensitively. A pattern without a
     * glob must be equal to the name; a pattern with one or more globs matches any run of
     * characters in their place (so {@code RE*AD} is legal in shape and simply matches nothing).
     */
    static boolean matches(String pattern, String name) {
        if( !pattern.contains(GLOB) )
            return pattern == name
        return toPattern(pattern).matcher(name).matches()
    }

    private static Pattern toPattern(String glob) {
        final sb = new StringBuilder()
        final String[] parts = glob.split(Pattern.quote(GLOB), -1)
        for( int i = 0; i < parts.length; i++ ) {
            if( i > 0 )
                sb.append('.*')
            if( parts[i] )
                sb.append(Pattern.quote(parts[i]))
        }
        return Pattern.compile(sb.toString())
    }
}
