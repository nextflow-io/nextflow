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
import groovy.transform.EqualsAndHashCode
import nextflow.exception.ScriptRuntimeException

/**
 * Expands the declared {@code tools} entries of an agent into the concrete set of tools the
 * agent may use.
 *
 * <p>The resolver is a <b>pure function</b> of (declared refs, available members): the members
 * of every family are supplied by the caller — the in-scope process names for
 * {@code nf:module_run}, the release-fixed leaf lists for {@code fs:} and {@code shell:} — so
 * the whole grammar can be exercised without a {@link nextflow.Session}, a script, or a runner.
 * A runner constraint (the {@code shell:} family needs a container boundary) likewise enters as
 * a <i>reason string</i> rather than as a runner name compared inside here.
 *
 * <p>Semantics, following the grammar:
 * <ul>
 *   <li><b>G3</b> — a ref names a node in a hierarchy: a leaf denotes one tool, a non-leaf its
 *       entire subtree, so {@code nf:module_run} is exactly {@code nf:module_run:*};</li>
 *   <li><b>G8</b> — selecting nothing is always an error, in five distinguishable flavours:
 *       malformed (raised by {@link ToolRef}), unknown family, missing explicit leaf, glob
 *       matching nothing, and non-leaf over an empty subtree. A directive is a declaration,
 *       not a filter;</li>
 *   <li><b>G9</b> — the result is an order-independent union. Overlapping refs are legal and
 *       idempotent, and the resolved list is emitted in <b>inventory</b> order (families in
 *       declaration-of-inventory order, members in the order the caller supplied them), never
 *       in the order the entries happen to appear in the directive.</li>
 * </ul>
 *
 * <p>The resolved tools carry the split the runner mapping needs: {@link ToolKind#BROKERED}
 * tools ({@code nf:module_run:X}) are executed by the driver as real Nextflow tasks and become
 * tool descriptors, whereas {@link ToolKind#NATIVE} tools ({@code fs:}, {@code shell:}) are
 * served by the runner itself and travel to it as bare names.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ToolRefResolver {

    // -- the standard inventory (§3). The names are the declaration-side ones; a leaf name is
    //    also its wire name (§4), which is why they are the bare Pi-baseline names.

    static final String NF_FAMILY = 'nf'

    static final String MODULE_RUN = 'module_run'

    static final String FS_FAMILY = 'fs'

    static final String SHELL_FAMILY = 'shell'

    /** The `fs:` leaves, in canonical order. Note there is no `exists`: `ls`/`read` cover it. */
    static final List<String> FS_TOOLS = List.of('read', 'write', 'edit', 'ls', 'grep', 'find')

    /** The `shell:` leaves. Its own family because it escapes the filesystem sandbox. */
    static final List<String> SHELL_TOOLS = List.of('bash')

    /** Who executes a resolved tool. */
    static enum ToolKind {
        /** Executed by the driver as a Nextflow task, wherever the agent runs. */
        BROKERED,
        /** Served by the runner's own tool of that name. */
        NATIVE
    }

    /**
     * A node of the tool hierarchy. A node with {@code null} members is a leaf, i.e. one tool;
     * a node with members is a group denoting its subtree. An <b>empty</b> member list is a
     * group with nothing in it — legal to model (a script with no processes in scope) and a
     * G8(e) error to select.
     */
    @CompileStatic
    static class ToolNode {
        final String name
        final List<ToolNode> members

        private ToolNode(String name, List<ToolNode> members) {
            this.name = name
            this.members = members
        }

        static ToolNode leaf(String name) { new ToolNode(name, null) }

        static ToolNode group(String name, List<ToolNode> members) { new ToolNode(name, members) }

        static List<ToolNode> leaves(Collection<String> names) {
            return names.collect { leaf(it) }
        }

        boolean isLeaf() { members == null }
    }

    /** A family: the root of one namespace, everything under it executed the same way. */
    @CompileStatic
    static class ToolFamily {
        final String name
        final ToolKind kind
        final List<ToolNode> members
        /** Why this family cannot be served here; {@code null} when it can. */
        final String unavailable
        /** Extra guidance appended to the G8(e) error when the family resolves to nothing. */
        final String emptyHint

        ToolFamily(String name, ToolKind kind, List<ToolNode> members, String unavailable = null, String emptyHint = null) {
            this.name = name
            this.kind = kind
            this.members = members
            this.unavailable = unavailable
            this.emptyHint = emptyHint
        }
    }

    /** One tool selected by the directive. */
    @CompileStatic
    @EqualsAndHashCode
    static class ResolvedTool {
        /** The canonical fully-qualified ref, e.g. {@code nf:module_run:GREET} or {@code fs:read}. */
        final String ref
        /** The wire name the model sees — the leaf name, never colon-bearing (§4). */
        final String name
        final ToolKind kind

        ResolvedTool(String ref, String name, ToolKind kind) {
            this.ref = ref
            this.name = name
            this.kind = kind
        }

        @Override
        String toString() { "${ref} -> ${name} (${kind})" }
    }

    /**
     * The resolved selection: every tool the directive selected, in inventory order, split by
     * who executes it. The two views are what the phases downstream consume — the brokered
     * names are wired into the tool bridge as descriptors, the native ones travel to the runner
     * beside them as bare names and never enter {@code toolSpecs}.
     */
    @CompileStatic
    static class Selection {
        private final List<ResolvedTool> tools

        Selection(List<ResolvedTool> tools) {
            this.tools = List.copyOf(tools)
        }

        List<ResolvedTool> getTools() { tools }

        boolean isEmpty() { tools.isEmpty() }

        private List<ResolvedTool> getBrokered() { tools.findAll { it.kind == ToolKind.BROKERED } }

        /** Named {@code runnerNative} because {@code native} is a reserved word. */
        List<ResolvedTool> getRunnerNative() { tools.findAll { it.kind == ToolKind.NATIVE } }

        /** The in-scope process names selected through {@code nf:module_run}. */
        List<String> getBrokeredNames() { getBrokered().collect { it.name } }

        /** The wire names of the runner-native tools, e.g. {@code [read, write, bash]}. */
        List<String> getNativeNames() { getRunnerNative().collect { it.name } }

        /**
         * The canonical refs of the runner-native tools, e.g. {@code [fs:read, shell:bash]}.
         * These are what the resume key folds in (a native tool has no descriptor to hash).
         */
        List<String> getNativeRefs() { getRunnerNative().collect { it.ref } }

        @Override
        String toString() { tools.toString() }
    }

    /** Per-ref bookkeeping, so a ref that selects nothing can say WHY it selected nothing. */
    private static class MatchStats {
        /** leaves contributed to the selection */
        int leaves
        /**
         * Whether the ref reached anything at all: a target node it matched, or a group along the
         * way that turned out to be empty. Both mean "the ref names something real that holds no
         * tools" — G8(e) — as opposed to naming nothing, so one flag serves where the two were
         * only ever read together.
         */
        boolean reachedNode
    }

    /** Family name -> family, iteration order = the order tools are emitted in. */
    private final Map<String, ToolFamily> families

    /** Prefixed to every error so the message names the agent; may be {@code null}. */
    private final String subject

    ToolRefResolver(String subject, Map<String, ToolFamily> families) {
        this.subject = subject
        this.families = families
    }

    /**
     * Build a resolver over the standard inventory: {@code nf:module_run} over the in-scope
     * processes, plus the release-fixed {@code fs:} and {@code shell:} leaves.
     *
     * @param subject          label prefixed to every error, e.g. {@code Agent `qc`}
     * @param processNames     the in-scope process names, i.e. the members of {@code nf:module_run}
     * @param shellUnavailable when non-{@code null}, the reason the {@code shell:} family cannot
     *                         be served by the selected runner; declaring any {@code shell:} ref
     *                         then fails with it. The runner is never named in here.
     */
    static ToolRefResolver standard(String subject, Collection<String> processNames, String shellUnavailable = null) {
        final families = new LinkedHashMap<String, ToolFamily>()
        families.put(NF_FAMILY, new ToolFamily(
                NF_FAMILY,
                ToolKind.BROKERED,
                List.of(ToolNode.group(MODULE_RUN, ToolNode.leaves(processNames ?: List.<String> of()))),
                null,
                'declare or `include` the processes the agent may run'))
        families.put(FS_FAMILY, new ToolFamily(FS_FAMILY, ToolKind.NATIVE, ToolNode.leaves(FS_TOOLS)))
        families.put(SHELL_FAMILY, new ToolFamily(SHELL_FAMILY, ToolKind.NATIVE, ToolNode.leaves(SHELL_TOOLS), shellUnavailable))
        return new ToolRefResolver(subject, families)
    }

    /**
     * Expand the declared entries into the resolved selection.
     *
     * @param declared the raw directive entries; each is resolved with {@code toString()} so an
     *                 interpolated entry is validated by its value (G1)
     * @return the resolved selection, empty only when nothing was declared (G7)
     * @throws ScriptRuntimeException on any malformed or zero-match entry (G8)
     */
    Selection resolve(Collection<?> declared) {
        // the canonical refs selected so far; a Set is what makes an overlapping ref idempotent
        // and the union order-independent (G9)
        final Set<String> selected = new HashSet<String>()
        if( declared ) {
            for( final entry : declared )
                select(entry?.toString(), selected)
        }
        // emit in INVENTORY order, not declaration order: the resolved set is a property of what
        // was selected, never of how it was written
        final List<ResolvedTool> result = new ArrayList<ResolvedTool>()
        for( final family : families.values() )
            emit(family.members, family.name, family.kind, selected, result)
        return new Selection(result)
    }

    /** Resolve one declared entry into {@code selected}, or throw explaining why it matched nothing. */
    private void select(String value, Set<String> selected) {
        ToolRef ref
        try {
            ref = ToolRef.parse(value)
        }
        catch( ScriptRuntimeException e ) {
            // G8(a) - re-thrown so the message names the agent as well as the ref
            throw fail(e.message)
        }
        final family = families.get(ref.family)
        // G8(b)
        if( family == null )
            throw fail("Tool `${ref}` belongs to the unknown family `${ref.family}` - known families are ${quoted(families.keySet())}")
        if( family.unavailable )
            throw fail("Tool `${ref}` is not available - ${family.unavailable}")

        final stats = new MatchStats()
        descend(family.members, ref, 1, family.name, selected, stats)
        if( stats.leaves > 0 )
            return
        // G8(c)/(d)/(e) - three distinct failures, told apart by how far the walk got.
        // An empty group anywhere along the way is reported as (e) whether the ref named it
        // directly or globbed below it, so `nf:module_run` and `nf:module_run:*` -- the same
        // selection by G3 -- also fail the same way.
        if( stats.reachedNode )
            throw fail("Tool `${ref}` selects nothing - it has no members${family.emptyHint ? '; ' + family.emptyHint : ''}")
        if( ref.globbed )
            throw fail("Tool pattern `${ref}` matches no tool - available: ${quoted(inventory(family))}; matching is case-sensitive")
        throw fail("Tool `${ref}` does not exist - available: ${quoted(inventory(family))}")
    }

    /**
     * Walk one level of the hierarchy matching {@code segments[idx]}. On reaching the terminal
     * segment the matched node is the target: a leaf contributes itself, a group its entire
     * subtree (G3). A ref deeper than the tree simply matches nothing.
     */
    private void descend(List<ToolNode> members, ToolRef ref, int idx, String prefix, Set<String> selected, MatchStats stats) {
        final segments = ref.segments
        final pattern = segments[idx]
        final terminal = idx == segments.size() - 1
        for( final node : members ) {
            if( !ToolRef.matches(pattern, node.name) )
                continue
            final path = "${prefix}:${node.name}".toString()
            if( terminal ) {
                stats.reachedNode = true
                collect(node, path, selected, stats)
            }
            else if( !node.isLeaf() ) {
                if( node.members.isEmpty() )
                    stats.reachedNode = true
                descend(node.members, ref, idx + 1, path, selected, stats)
            }
        }
    }

    /** Add every leaf of the target node's subtree; the node itself when it is a leaf. */
    private void collect(ToolNode node, String path, Set<String> selected, MatchStats stats) {
        if( node.isLeaf() ) {
            selected.add(path)
            stats.leaves++
            return
        }
        for( final child : node.members )
            collect(child, "${path}:${child.name}".toString(), selected, stats)
    }

    /** Emit the selected leaves of a subtree, preserving the order the members were declared in. */
    private static void emit(List<ToolNode> members, String prefix, ToolKind kind, Set<String> selected, List<ResolvedTool> out) {
        for( final node : members ) {
            final path = "${prefix}:${node.name}".toString()
            if( node.isLeaf() ) {
                if( selected.contains(path) )
                    out.add(new ResolvedTool(path, node.name, kind))
            }
            else {
                emit(node.members, path, kind, selected, out)
            }
        }
    }

    /** Every selectable ref of a family, for the "available:" part of a zero-match error. */
    private static List<String> inventory(ToolFamily family) {
        final out = new LinkedHashSet<String>()
        collectAll(family.members, family.name, out)
        return new ArrayList<String>(out)
    }

    private static void collectAll(List<ToolNode> members, String prefix, Set<String> out) {
        for( final node : members ) {
            final path = "${prefix}:${node.name}".toString()
            if( node.isLeaf() )
                out.add(path)
            else
                collectAll(node.members, path, out)
        }
    }

    private static String quoted(Collection<String> values) {
        return values ? values.collect { "`${it}`" }.join(', ') : '(none)'
    }

    private ScriptRuntimeException fail(String message) {
        return new ScriptRuntimeException(subject ? "${subject}: ${message}".toString() : message)
    }
}
