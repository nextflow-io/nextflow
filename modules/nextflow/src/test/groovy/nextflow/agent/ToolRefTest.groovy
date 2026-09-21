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

import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Grammar-level tests: the SHAPE of a single tool ref, independent of which families or tools
 * happen to exist. Whether a well-formed ref actually selects anything is {@link ToolRefResolver}'s
 * business and is covered by {@code ToolRefResolverTest}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ToolRefTest extends Specification {

    @Unroll
    def 'should accept the well-formed ref #REF'() {
        when:
        def ref = ToolRef.parse(REF)
        then:
        ref.ref == REF
        ref.segments == SEGMENTS
        ref.family == SEGMENTS[0]
        ref.globbed == GLOBBED

        where:
        REF                             | SEGMENTS                              | GLOBBED
        'nf:module_run'                 | ['nf', 'module_run']                  | false
        'nf:module_run:SAMTOOLS_SORT'   | ['nf', 'module_run', 'SAMTOOLS_SORT'] | false
        'nf:module_run:SAMTOOLS_*'      | ['nf', 'module_run', 'SAMTOOLS_*']    | true
        'nf:module_run:*'               | ['nf', 'module_run', '*']             | true
        'nf:*'                          | ['nf', '*']                           | true
        'fs:*'                          | ['fs', '*']                           | true
        'fs:read'                       | ['fs', 'read']                        | false
        'shell:bash'                    | ['shell', 'bash']                     | false
        // shape-legal even though it can only ever match nothing (that is a G8(d) failure,
        // raised by the resolver, not a syntax error)
        'fs:RE*AD'                      | ['fs', 'RE*AD']                       | true
        // a family that does not exist is still a well-formed ref; G8(b) rejects it later
        'mcp:github:get_issue'          | ['mcp', 'github', 'get_issue']        | false
    }

    @Unroll
    def 'should reject #REF because #WHY'() {
        when:
        ToolRef.parse(REF)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Invalid tool reference `${REF}`")
        e.message.contains(REASON)

        where:
        REF                       | WHY                                      | REASON
        // -- G1: no fallthrough, the four legacy non-capability forms are all errors
        'greet'                   | 'a bare process name is not namespaced'  | 'must be namespaced'
        './mod.nf'                | 'a local module path is not namespaced'  | 'must be namespaced'
        '/abs/path/mod.nf'        | 'an absolute path is not namespaced'     | 'must be namespaced'
        'nf-core/skesa'           | 'a registry ref is not namespaced'       | 'must be namespaced'
        'module_run'              | 'the old capability string is bare'      | 'must be namespaced'
        'filesystem'              | 'the old capability string is bare'      | 'must be namespaced'
        // -- G5: a glob with no family means "every tool that exists"
        '*'                       | 'the glob is unanchored'                 | 'must be anchored to a tool family'
        'SAMTOOLS_*'              | 'the glob is unanchored'                 | 'must be anchored to a tool family'
        // -- G4: globs are trailing only
        'nf*:module_run'          | 'the family carries a glob'              | '`nf*` cannot contain a glob'
        'nf:*:SAMTOOLS_SORT'      | 'an intermediate segment carries a glob' | '`*` cannot contain a glob'
        // -- G6: there is no exclude operator
        '!fs:bash'                | 'negation is not supported'              | 'exclusions are not supported'
        '!nf:module_run'          | 'negation is not supported'              | 'exclusions are not supported'
        // -- G2: at least two segments, none empty, never trimmed
        'nf::read'                | 'the middle segment is empty'            | 'segment 2 is empty'
        ':read'                   | 'the family is empty'                    | 'segment 1 is empty'
        'fs:'                     | 'the last segment is empty'              | 'segment 2 is empty'
        'nf'                      | 'a single segment is not a ref'          | 'must be namespaced'
        'fs:read '                | 'the entry is never trimmed'             | 'is not a legal segment'
        ' fs:read'                | 'the entry is never trimmed'             | 'is not a legal segment'
        'fs:re ad'                | 'a segment cannot contain a space'       | 'is not a legal segment'
        '~/Projects/nf:dev'       | 'a path is not a legal family'           | 'is not a legal segment'
        'nf:module.run'           | 'a dot is not a legal segment character' | 'is not a legal segment'
    }

    def 'should reject an empty entry'() {
        when:
        ToolRef.parse(VALUE)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('cannot be empty')

        where:
        VALUE << [null, '']
    }

    @Unroll
    def 'should match #PATTERN against #NAME case-sensitively'() {
        expect:
        ToolRef.matches(PATTERN, NAME) == EXPECTED

        where:
        PATTERN         | NAME             | EXPECTED
        'read'          | 'read'           | true
        'read'          | 'Read'           | false
        '*'             | 'read'           | true
        '*'             | ''               | true
        'SAMTOOLS_*'    | 'SAMTOOLS_SORT'  | true
        'SAMTOOLS_*'    | 'SAMTOOLS'       | false
        'SAMTOOLS_*'    | 'BWA_MEM'        | false
        // G2: matching is case-sensitive in every segment, globs included
        'samtools_*'    | 'SAMTOOLS_SORT'  | false
        'RE*AD'         | 'read'           | false
        'RE*AD'         | 'REAAD'          | true
        'RE*AD'         | 'REALLY_BAD'     | true
        // the glob is not a substring search: the pattern must span the WHOLE name
        'RE*AD'         | 'REALLY_BADLY'   | false
        'RE*AD'         | 'THREAD'         | false
        're*d'          | 'read'           | true
        '*_SORT'        | 'SAMTOOLS_SORT'  | true
        // `-` is a legal segment character and must not be read as a regex range
        'a-*'           | 'a-b'            | true
    }
}
