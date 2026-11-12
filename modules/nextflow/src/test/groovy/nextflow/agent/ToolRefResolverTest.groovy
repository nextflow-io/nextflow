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

import nextflow.agent.ToolRefResolver.ToolKind
import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Resolution tests: what a well-formed ref actually selects, given the members available.
 * The resolver is a pure function of (declared refs, available members), so every case below
 * runs without a session, a script or a runner.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ToolRefResolverTest extends Specification {

    static final List<String> PROCESSES = ['GREET', 'SAMTOOLS_SORT', 'SAMTOOLS_INDEX']

    private static ToolRefResolver resolver(List<String> processes = PROCESSES, String shellUnavailable = null) {
        return ToolRefResolver.standard('Agent `qc`', processes, shellUnavailable)
    }

    // ------------------------------------------------------------------ acceptance

    @Unroll
    def 'should expand #DECLARED to #NAMES'() {
        when:
        def selection = resolver().resolve(DECLARED)
        then:
        selection.tools.collect { it.name } == NAMES
        selection.tools.collect { it.ref } == REFS

        where:
        DECLARED                             | NAMES                                          | REFS
        // -- G3: a non-leaf denotes its entire subtree, so `nf:module_run` is `nf:module_run:*`
        ['nf:module_run']                    | ['GREET', 'SAMTOOLS_SORT', 'SAMTOOLS_INDEX']   | ['nf:module_run:GREET', 'nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
        ['nf:module_run:*']                  | ['GREET', 'SAMTOOLS_SORT', 'SAMTOOLS_INDEX']   | ['nf:module_run:GREET', 'nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
        // today `nf:*` and `nf:module_run` coincide because module_run is the family's only
        // member -- an inventory coincidence, not an alias
        ['nf:*']                             | ['GREET', 'SAMTOOLS_SORT', 'SAMTOOLS_INDEX']   | ['nf:module_run:GREET', 'nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
        ['nf:module_run:GREET']              | ['GREET']                                      | ['nf:module_run:GREET']
        // -- G4: a trailing glob, matched case-sensitively
        ['nf:module_run:SAMTOOLS_*']         | ['SAMTOOLS_SORT', 'SAMTOOLS_INDEX']            | ['nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
        // -- G5: the family root anchors a glob over release-fixed membership
        ['fs:*']                             | ['read', 'write', 'edit', 'ls', 'grep', 'find'] | ['fs:read', 'fs:write', 'fs:edit', 'fs:ls', 'fs:grep', 'fs:find']
        ['fs:read']                          | ['read']                                       | ['fs:read']
        ['fs:read', 'fs:write']              | ['read', 'write']                              | ['fs:read', 'fs:write']
        ['shell:bash']                       | ['bash']                                       | ['shell:bash']
        ['shell:*']                          | ['bash']                                       | ['shell:bash']
        // -- a mixed directive: brokered first because `nf` is the first family of the inventory
        ['fs:*', 'nf:module_run:GREET']      | ['GREET', 'read', 'write', 'edit', 'ls', 'grep', 'find'] | ['nf:module_run:GREET', 'fs:read', 'fs:write', 'fs:edit', 'fs:ls', 'fs:grep', 'fs:find']
    }

    def 'should carry the runner split on each resolved tool'() {
        when:
        def selection = resolver().resolve(['nf:module_run:GREET', 'fs:read', 'shell:bash'])
        then:
        selection.tools.collect { it.kind } == [ToolKind.BROKERED, ToolKind.NATIVE, ToolKind.NATIVE]
        and: 'the brokered half is the process names the driver wires as tasks'
        selection.brokeredNames == ['GREET']
        and: 'the native half is bare wire names for the runner, plus the refs the cache key needs'
        selection.nativeNames == ['read', 'bash']
        selection.nativeRefs == ['fs:read', 'shell:bash']
    }

    def 'should resolve an interpolated entry by its value'() {
        given:
        def proc = 'GREET'
        when:
        def selection = resolver().resolve(["nf:module_run:${proc}"])
        then:
        selection.tools.collect { it.name } == ['GREET']
    }

    def 'should return an empty selection when nothing is declared'() {
        expect:
        resolver().resolve(DECLARED).isEmpty()
        where:
        DECLARED << [null, [], [] as Set]
    }

    // ------------------------------------------------------------------ G9 union

    @Unroll
    def 'should collapse the overlapping refs #DECLARED'() {
        when:
        def selection = resolver().resolve(DECLARED)
        then:
        selection.tools.collect { it.ref } == REFS

        where:
        DECLARED                                        | REFS
        ['fs:*', 'fs:read']                             | ['fs:read', 'fs:write', 'fs:edit', 'fs:ls', 'fs:grep', 'fs:find']
        ['fs:read', 'fs:read']                          | ['fs:read']
        ['nf:module_run', 'nf:module_run:GREET']        | ['nf:module_run:GREET', 'nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
        ['nf:*', 'nf:module_run:SAMTOOLS_*']            | ['nf:module_run:GREET', 'nf:module_run:SAMTOOLS_SORT', 'nf:module_run:SAMTOOLS_INDEX']
    }

    @Unroll
    def 'should resolve the same set whatever the entry order'() {
        expect:
        resolver().resolve(A).tools == resolver().resolve(B).tools

        where:
        A                                          | B
        ['fs:read', 'fs:write']                    | ['fs:write', 'fs:read']
        ['fs:*', 'nf:module_run']                  | ['nf:module_run', 'fs:*']
        ['nf:module_run:GREET', 'fs:*', 'shell:bash'] | ['shell:bash', 'fs:*', 'nf:module_run:GREET']
        ['fs:*', 'fs:read']                        | ['fs:read', 'fs:*']
    }

    // ------------------------------------------------------------------ G8 zero-match

    def 'should reject a malformed ref naming the agent'() {
        when: 'G8(a)'
        resolver().resolve(['greet'])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.startsWith('Agent `qc`: ')
        e.message.contains('Invalid tool reference `greet`')
        e.message.contains('must be namespaced')
    }

    def 'should reject an unknown family'() {
        when: 'G8(b)'
        resolver().resolve(['mcp:github:get_issue'])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('unknown family `mcp`')
        e.message.contains('known families are `nf`, `fs`, `shell`')
    }

    def 'should reject an explicit leaf that does not exist'() {
        when: 'G8(c)'
        resolver().resolve(['nf:module_run:MISSING'])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Tool `nf:module_run:MISSING` does not exist')
        e.message.contains('`nf:module_run:GREET`')
    }

    def 'should reject an explicit leaf whose case does not match'() {
        when: 'G2 matching is case-sensitive, so this is a G8(c)'
        resolver().resolve(['nf:module_run:greet'])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Tool `nf:module_run:greet` does not exist')
    }

    @Unroll
    def 'should reject the glob #REF that matches nothing'() {
        when: 'G8(d)'
        resolver().resolve([REF])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Tool pattern `${REF}` matches no tool")
        e.message.contains('case-sensitive')

        where:
        REF << [
            'nf:module_run:ZZZ_*',
            // shape-legal, but matches none of the six fs leaves
            'fs:RE*AD',
            // the process segment keeps process-name case
            'nf:module_run:samtools_*' ]
    }

    @Unroll
    def 'should reject #REF over an empty subtree'() {
        given: 'a script with no processes in scope'
        def resolver = resolver([])
        when: 'G8(e)'
        resolver.resolve([REF])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Tool `${REF}` selects nothing")
        e.message.contains('include')

        where:
        REF << ['nf:module_run', 'nf:*', 'nf:module_run:*']
    }

    def 'should keep the five zero-match failures distinguishable'() {
        given:
        def r = resolver()
        when:
        r.resolve([BAD])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains(EXPECTED)

        where:
        BAD                     | EXPECTED
        'greet'                 | 'Invalid tool reference'
        'mcp:x'                 | 'unknown family'
        'nf:module_run:MISSING' | 'does not exist'
        'nf:module_run:ZZZ_*'   | 'matches no tool'
        'fs:read:deeper'        | 'does not exist'
    }

    // ------------------------------------------------------------------ runner constraint

    def 'should reject the shell family when the runner cannot serve it'() {
        given:
        def r = resolver(PROCESSES, 'the `pi` runner is required')
        when:
        r.resolve([REF])
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Tool `${REF}` is not available")
        e.message.contains('the `pi` runner is required')

        where:
        REF << ['shell:bash', 'shell:*']
    }

    def 'should still serve the other families when shell is unavailable'() {
        when:
        def selection = resolver(PROCESSES, 'nope').resolve(['fs:read', 'nf:module_run:GREET'])
        then:
        selection.tools.collect { it.ref } == ['nf:module_run:GREET', 'fs:read']
    }
}
