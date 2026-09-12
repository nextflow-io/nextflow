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

import java.nio.file.Files
import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Unit tests for the six {@code fs:} tools — {@code read}, {@code write}, {@code edit},
 * {@code ls}, {@code grep}, {@code find} — as served in the driver JVM by
 * {@link ModuleToolBridge}, plus their {@link FilesystemTools} descriptors.
 *
 * <p>The bridge is built with NO modules, so every call exercises the filesystem dispatch
 * path alone. Each tool is covered by a happy path and by both sandbox-escape shapes: a
 * {@code ..} traversal out of the work dir, and a symlink INSIDE the work dir whose target
 * is outside it — the case a lexical path check would miss.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilesystemToolsTest extends Specification {

    @TempDir Path workDir

    /** A directory outside the sandbox, created lazily by {@link #outsideFile}. */
    private Path outside

    def cleanup() {
        ModuleToolBridge.clearContext()
    }

    /** A bridge with no modules serving the given fs: tools (all six by default). */
    private ModuleToolBridge bridge(Collection<String> fsTools = FilesystemTools.NAMES) {
        new ModuleToolBridge(
            Collections.<WiredModuleTool>emptyList(),
            fsTools )
    }

    /** Invoke a tool and parse its JSON result. */
    private Map call(ModuleToolBridge bridge, String tool, Map args) {
        return (Map) new JsonSlurper().parseText(bridge.call(tool, JsonOutput.toJson(args)))
    }

    /** Create a file OUTSIDE the sandbox (a sibling of the work dir) and return it. */
    private Path outsideFile(String name, String content) {
        if( outside == null )
            outside = Files.createDirectories(workDir.getParent().resolve('fs-outside-' + System.nanoTime()))
        return Files.write(outside.resolve(name), content.getBytes('UTF-8'))
    }

    /** A symlink inside the work dir pointing at a target outside it. */
    private Path escapingLink(String linkName, Path target) {
        return Files.createSymbolicLink(workDir.resolve(linkName), target)
    }

    private void sandbox() {
        ModuleToolBridge.setContext(new DispatchContext(workDir))
    }

    // =========================================================================
    // descriptors (§4 wire names)
    // =========================================================================

    def 'should advertise the six fs tools under their bare wire names'() {
        when:
        def names = bridge().filesystemDescriptors()*.name
        then:
        names == ['read','write','edit','ls','grep','find']
        and: 'the legacy aggregate tool is gone for good'
        !names.contains('filesystem')
        and: 'no wire name is colon-bearing, and all are OpenAI-legal'
        names.every { it ==~ /[a-zA-Z0-9_-]{1,64}/ }
        and: '§5: a runner-native tool is never a brokered descriptor'
        bridge().descriptors().isEmpty()
    }

    def 'should advertise only the selected subset, in canonical order'() {
        expect: 'the descriptor order follows the family inventory, never the declaration order'
        bridge(['grep','read']).filesystemDescriptors()*.name == ['read','grep']
        bridge([]).filesystemDescriptors().isEmpty()
        and: 'the sandbox flag follows the selection'
        bridge(['read']).filesystemEnabled
        !bridge([]).filesystemEnabled
    }

    def 'should give every grammar leaf of the fs family a descriptor'() {
        expect: 'the selectable set and the servable set cannot drift apart'
        FilesystemTools.NAMES == ToolRefResolver.FS_TOOLS
        ToolRefResolver.FS_TOOLS.every { FilesystemTools.descriptor(it) != null }
    }

    def 'should reject a descriptor request for an unknown tool'() {
        when:
        FilesystemTools.descriptor('exists')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('exists')
    }

    def 'should state the sandbox boundary in every description'() {
        expect:
        FilesystemTools.descriptors(FilesystemTools.NAMES).every {
            it.description.toLowerCase().contains('sandbox') || it.description.contains('work dir')
        }
    }

    def 'should not route a tool name that was not selected'() {
        given: 'an agent that declared fs:read only'
        def bridge = bridge(['read'])
        sandbox()
        when: 'the model calls a filesystem tool it was never given'
        def result = call(bridge, 'grep', [pattern: 'x'])
        then: 'it is an unknown tool, not a hijacked filesystem call'
        result.error.contains('Unknown agent tool `grep`')
    }

    // =========================================================================
    // no dispatch context
    // =========================================================================

    def 'should fail every fs tool without a sandbox context'() {
        given:
        def bridge = bridge()
        expect:
        call(bridge, tool, args).error.contains('no sandbox context')
        where:
        tool    | args
        'read'  | [path: 'a.txt']
        'write' | [path: 'a.txt', content: 'x']
        'edit'  | [path: 'a.txt', old_string: 'a', new_string: 'b']
        'ls'    | [:]
        'grep'  | [pattern: 'x']
        'find'  | [pattern: '*.txt']
    }

    // =========================================================================
    // read
    // =========================================================================

    def 'read should return the content of a text file'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('data.txt'), 'the content'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'read', [path: 'data.txt'])
        then:
        result.content == 'the content'
    }

    def 'read should return a path handle for a non text-like file'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('reads.fastq'), 'ACGT'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'read', [path: 'reads.fastq'])
        then: 'the opaque-path contract: chainable, but the bytes are not inlined - and never mislabelled as content'
        result.content == null
        result.error == null
        result.path.endsWith('reads.fastq')
        result.note.contains('not inlined')
    }

    def 'read should not label a binary file as content'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('sneaky.txt'), ['a'.bytes, [0 as byte] as byte[], 'b'.bytes].flatten() as byte[])
        sandbox()
        when:
        def result = call(bridge, 'read', [path: 'sneaky.txt'])
        then:
        result.content == null
        result.note.contains('binary')
    }

    def 'read should honour the inline size cap'() {
        given:
        def bridge = bridge()
        bridge.setMaxInlineBytes(8)
        Files.write(workDir.resolve('big.txt'), ('x' * 100).bytes)
        sandbox()
        when:
        def result = call(bridge, 'read', [path: 'big.txt'])
        then:
        result.note.contains('content not inlined')
    }

    def 'read should reject a missing file and a directory'() {
        given:
        def bridge = bridge()
        Files.createDirectories(workDir.resolve('sub'))
        sandbox()
        expect:
        call(bridge, 'read', [path: 'nope.txt']).error.contains('file not found')
        call(bridge, 'read', [path: 'sub']).error.contains('is a directory')
    }

    def 'read should refuse a .. traversal out of the sandbox'() {
        given:
        def bridge = bridge()
        outsideFile('secret.txt', 'topsecret')
        sandbox()
        when:
        def result = call(bridge, 'read', [path: '../secret.txt'])
        then:
        result.error.contains('outside sandbox')
    }

    def 'read should refuse a symlink whose target is outside the sandbox'() {
        given:
        def bridge = bridge()
        escapingLink('link.txt', outsideFile('secret.txt', 'topsecret'))
        sandbox()
        when:
        def result = call(bridge, 'read', [path: 'link.txt'])
        then:
        result.error.contains('outside sandbox')
        result.content == null
    }

    // =========================================================================
    // write
    // =========================================================================

    def 'write should create a file and its parent dirs'() {
        given:
        def bridge = bridge()
        sandbox()
        when:
        def result = call(bridge, 'write', [path: 'nested/out.txt', content: 'hello world'])
        then:
        result.bytes == 11
        workDir.resolve('nested/out.txt').text == 'hello world'
    }

    def 'write should refuse a .. traversal out of the sandbox'() {
        given:
        def bridge = bridge()
        sandbox()
        when:
        def result = call(bridge, 'write', [path: '../evil.txt', content: 'bad'])
        then:
        result.error.contains('outside sandbox')
        !Files.exists(workDir.getParent().resolve('evil.txt'))
    }

    def 'write should refuse a symlink whose target is outside the sandbox'() {
        given:
        def bridge = bridge()
        def target = outsideFile('victim.txt', 'original')
        escapingLink('victim.txt', target)
        sandbox()
        when:
        def result = call(bridge, 'write', [path: 'victim.txt', content: 'overwritten'])
        then:
        result.error.contains('outside sandbox')
        target.text == 'original'
    }

    def 'write should refuse a whitelisted read-only path'() {
        given: 'a module output outside the work dir, readable but never writable'
        def bridge = bridge()
        def output = outsideFile('module-out.txt', 'produced')
        def ctx = new DispatchContext(workDir)
        ctx.addReadablePath(output)
        ModuleToolBridge.setContext(ctx)
        expect:
        call(bridge, 'read', [path: output.toAbsolutePath().toString()]).content == 'produced'
        call(bridge, 'write', [path: output.toAbsolutePath().toString(), content: 'x']).error.contains('outside sandbox')
        output.text == 'produced'
    }

    def 'write should refuse a missing content instead of truncating the file'() {
        given: 'an existing file the model must not destroy by omitting an argument'
        def bridge = bridge()
        Files.write(workDir.resolve('keep.txt'), 'precious'.bytes)
        sandbox()

        expect: 'a missing (or null) content is an error, exactly as edit treats new_string'
        call(bridge, 'write', [path: 'keep.txt']).error.contains('missing required argument: content')
        call(bridge, 'write', [path: 'keep.txt', content: null]).error.contains('missing required argument: content')
        workDir.resolve('keep.txt').text == 'precious'

        and: 'an empty file is still writable -- it just has to be asked for'
        call(bridge, 'write', [path: 'keep.txt', content: '']).bytes == 0
        workDir.resolve('keep.txt').text == ''
    }

    /**
     * The module-output whitelist is what makes {@code read} usable at all across a
     * {@code nf:module_run} call: the produced FILE becomes readable, its siblings do not.
     * Covered here because it is the {@code fs:} side of that contract.
     */
    def 'read should honour the module-output whitelist, file by file'() {
        given:
        def bridge = bridge()
        def output = outsideFile('result.txt', 'module output content')
        def sibling = outsideFile('sibling-secret.txt', 'not an output')
        def ctx = new DispatchContext(workDir)
        ModuleToolBridge.setContext(ctx)
        when: 'a successful module call reports the produced path'
        ModuleToolBridge.whitelistOutputDirs([output: output.toAbsolutePath().toString()])
        then: 'the file is whitelisted, its parent dir is not'
        ctx.readablePaths.contains(output)
        !ctx.readablePaths.contains(output.getParent())
        and:
        call(bridge, 'read', [path: output.toAbsolutePath().toString()]).content == 'module output content'
        call(bridge, 'read', [path: sibling.toAbsolutePath().toString()]).error.contains('outside sandbox')
    }

    // =========================================================================
    // edit
    // =========================================================================

    def 'edit should replace a unique occurrence'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('cfg.txt'), 'alpha\nbeta\ngamma\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'edit', [path: 'cfg.txt', old_string: 'beta', new_string: 'BETA'])
        then:
        result.replacements == 1
        workDir.resolve('cfg.txt').text == 'alpha\nBETA\ngamma\n'
    }

    def 'edit should delete the matched text when new_string is empty'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('cfg.txt'), 'keep\ndrop\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'edit', [path: 'cfg.txt', old_string: 'drop\n', new_string: ''])
        then:
        result.replacements == 1
        workDir.resolve('cfg.txt').text == 'keep\n'
    }

    def 'edit should refuse a non-unique old_string instead of editing the first match'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('dup.txt'), 'x = 1\ny = 2\nx = 1\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'edit', [path: 'dup.txt', old_string: 'x = 1', new_string: 'x = 9'])
        then: 'the count is reported and BOTH ways out are named'
        result.error.contains('occurs 2 times')
        result.error.contains('replace_all')
        and: 'the file is untouched - no silent first-match edit'
        workDir.resolve('dup.txt').text == 'x = 1\ny = 2\nx = 1\n'
    }

    def 'edit should replace every occurrence with replace_all'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('dup.txt'), 'x = 1\ny = 2\nx = 1\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'edit', [path: 'dup.txt', old_string: 'x = 1', new_string: 'x = 9', replace_all: true])
        then:
        result.replacements == 2
        workDir.resolve('dup.txt').text == 'x = 9\ny = 2\nx = 9\n'
    }

    def 'edit should report a missing match, an identical pair and a missing argument'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('cfg.txt'), 'alpha\n'.bytes)
        sandbox()
        expect:
        call(bridge, 'edit', [path: 'cfg.txt', old_string: 'zeta', new_string: 'x']).error.contains('no match for old_string')
        call(bridge, 'edit', [path: 'cfg.txt', old_string: 'alpha', new_string: 'alpha']).error.contains('identical')
        call(bridge, 'edit', [path: 'cfg.txt', old_string: 'alpha']).error.contains('new_string')
        call(bridge, 'edit', [path: 'cfg.txt', new_string: 'x']).error.contains('old_string')
        and: 'a literal old_string is never read as a regular expression'
        call(bridge, 'edit', [path: 'cfg.txt', old_string: 'a.p', new_string: 'x']).error.contains('no match for old_string')
    }

    def 'edit should refuse a .. traversal and a symlink out of the sandbox'() {
        given:
        def bridge = bridge()
        def target = outsideFile('victim.txt', 'original')
        escapingLink('victim.txt', target)
        sandbox()
        expect:
        call(bridge, 'edit', [path: '../victim.txt', old_string: 'original', new_string: 'hacked']).error.contains('outside sandbox')
        call(bridge, 'edit', [path: 'victim.txt', old_string: 'original', new_string: 'hacked']).error.contains('outside sandbox')
        target.text == 'original'
    }

    // =========================================================================
    // ls
    // =========================================================================

    def 'ls should list the immediate entries with their type and size'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'abc'.bytes)
        Files.createDirectories(workDir.resolve('sub'))
        Files.write(workDir.resolve('sub/deep.txt'), 'nested'.bytes)
        sandbox()
        when: 'the path argument is omitted it defaults to the work dir'
        def result = call(bridge, 'ls', [:])
        then:
        result.entries.find { it.name=='a.txt' } == [name: 'a.txt', type: 'file', size: 3]
        result.entries.find { it.name=='sub' } == [name: 'sub', type: 'dir']
        and: 'it is not recursive'
        !result.entries.any { it.name=='deep.txt' }
    }

    def 'ls should reject a missing dir and a file'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'abc'.bytes)
        sandbox()
        expect:
        call(bridge, 'ls', [path: 'nope']).error.contains('directory not found')
        call(bridge, 'ls', [path: 'a.txt']).error.contains('not a directory')
    }

    def 'ls should refuse a .. traversal and a symlink out of the sandbox'() {
        given:
        def bridge = bridge()
        outsideFile('secret.txt', 'topsecret')
        Files.createSymbolicLink(workDir.resolve('linkdir'), outside)
        sandbox()
        expect:
        call(bridge, 'ls', [path: '..']).error.contains('outside sandbox')
        call(bridge, 'ls', [path: 'linkdir']).error.contains('outside sandbox')
    }

    def 'ls should not report the type or size of an entry pointing outside the sandbox'() {
        given: 'a symlink INSIDE the work dir whose target is a file outside it'
        def bridge = bridge()
        def secret = outsideFile('secret.txt', 'a very specific number of bytes')
        escapingLink('innocent.txt', secret)
        Files.write(workDir.resolve('own.txt'), 'abc'.bytes)
        sandbox()

        when: 'the PARENT is listed -- `ls` on the link itself is already refused'
        def result = call(bridge, 'ls', [:])
        def escaping = result.entries.find { it.name == 'innocent.txt' }

        then: 'the entry is reported, so the later read refusal is not inexplicable'
        escaping != null
        and: 'but neither its kind nor its exact size leaks, both of which stat THROUGH the link'
        escaping == [name: 'innocent.txt', type: 'link']
        and: 'an entry genuinely inside the sandbox is unaffected'
        result.entries.find { it.name == 'own.txt' } == [name: 'own.txt', type: 'file', size: 3]
    }

    // =========================================================================
    // grep
    // =========================================================================

    def 'grep should find matching lines recursively with file and line number'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'nothing\nhello world\n'.bytes)
        Files.createDirectories(workDir.resolve('sub'))
        Files.write(workDir.resolve('sub/b.txt'), 'hello again\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'grep', [pattern: 'hello'])
        then:
        result.count == 2
        result.truncated == false
        result.matches.find { it.text=='hello world' }.line == 2
        result.matches.find { it.text=='hello again' }.file.endsWith('sub/b.txt')
    }

    def 'grep should accept a regular expression and a case_insensitive flag'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'Sample_01\nsample_02\nother\n'.bytes)
        sandbox()
        expect:
        call(bridge, 'grep', [pattern: '^sample_\\d+$']).count == 1
        call(bridge, 'grep', [pattern: '^sample_\\d+$', case_insensitive: true]).count == 2
        call(bridge, 'grep', [pattern: '[unclosed']).error.contains('invalid regular expression')
    }

    def 'grep should restrict the search with the include glob'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'needle\n'.bytes)
        Files.write(workDir.resolve('b.log'), 'needle\n'.bytes)
        sandbox()
        expect:
        call(bridge, 'grep', [pattern: 'needle']).count == 2
        call(bridge, 'grep', [pattern: 'needle', include: '*.log']).count == 1
    }

    def 'grep should cap the results and say so'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('many.txt'), (1..20).collect { "match ${it}" }.join('\n').bytes)
        sandbox()
        when:
        def result = call(bridge, 'grep', [pattern: 'match', max_results: 3])
        then: 'the cap travels back so the model knows it did not see everything'
        result.count == 3
        result.truncated == true
        result.limit == 3
        result.truncated_reason == FilesystemTools.TRUNCATED_RESULTS
        when: 'the cap is left to the default'
        def full = call(bridge, 'grep', [pattern: 'match'])
        then:
        full.count == 20
        full.truncated == false
        full.limit == FilesystemTools.DEFAULT_MAX_RESULTS
        and: 'a complete search never carries a reason to explain away'
        full.truncated_reason == null
    }

    /**
     * {@code max_results} bounds the ANSWER; this bounds the WORK. The two are reported
     * differently on purpose: "you saw the first N matches" and "the search stopped early" are
     * opposite instructions — narrow the pattern versus narrow the root — and a model that cannot
     * tell them apart reads an incomplete search as an exhaustive one that found nothing.
     */
    def 'grep should bound the search itself and name that as the reason'() {
        given: 'more candidate files than one call may read'
        def bridge = bridge()
        final over = FilesystemTools.MAX_GREP_FILES + 1
        for( int i=0; i<over; i++ )
            Files.write(workDir.resolve("f${i}.txt"), 'needle\n'.bytes)
        sandbox()

        when: 'a pattern that matches NOTHING -- the case the result limit can never bound'
        def result = call(bridge, 'grep', [pattern: 'no-such-token'])

        then: 'the file budget stops it, and says so distinctly'
        result.count == 0
        result.truncated == true
        result.truncated_reason == FilesystemTools.TRUNCATED_SEARCH
        result.files_scanned == FilesystemTools.MAX_GREP_FILES
    }

    def 'grep should clamp an over-large cap to the hard maximum'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'match\n'.bytes)
        sandbox()
        expect:
        call(bridge, 'grep', [pattern: 'match', max_results: 999999]).limit == FilesystemTools.MAX_MAX_RESULTS
        call(bridge, 'grep', [pattern: 'match', max_depth: 999999]).max_depth == FilesystemTools.MAX_MAX_DEPTH
    }

    def 'grep should cap the walk depth'() {
        given:
        def bridge = bridge()
        Files.createDirectories(workDir.resolve('l1/l2'))
        Files.write(workDir.resolve('top.txt'), 'needle\n'.bytes)
        Files.write(workDir.resolve('l1/l2/deep.txt'), 'needle\n'.bytes)
        sandbox()
        expect:
        call(bridge, 'grep', [pattern: 'needle', max_depth: 1]).count == 1
        call(bridge, 'grep', [pattern: 'needle']).count == 2
    }

    def 'grep should skip a binary file'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('bin.dat'), ['needle'.bytes, [0 as byte] as byte[]].flatten() as byte[])
        Files.write(workDir.resolve('ok.txt'), 'needle\n'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'grep', [pattern: 'needle'])
        then:
        result.count == 1
        result.files_scanned == 1
    }

    def 'grep should refuse a root outside the sandbox and skip an escaping symlink'() {
        given:
        def bridge = bridge()
        escapingLink('link.txt', outsideFile('secret.txt', 'topsecret'))
        Files.write(workDir.resolve('inside.txt'), 'topsecret is also here\n'.bytes)
        sandbox()
        when: 'the search root escapes the sandbox'
        def escaped = call(bridge, 'grep', [pattern: 'topsecret', path: '..'])
        then:
        escaped.error.contains('outside sandbox')
        when: 'a symlinked file inside the work dir points out of it'
        def result = call(bridge, 'grep', [pattern: 'topsecret'])
        then: 'only the genuinely-inside file is searched'
        result.count == 1
        result.matches[0].file.endsWith('inside.txt')
        !result.matches.any { it.file.contains('link.txt') }

        when: 'the very same target is whitelisted as a module output'
        def ctx = new DispatchContext(workDir)
        ctx.addReadablePath(outside.resolve('secret.txt'))
        ModuleToolBridge.setContext(ctx)
        def allowed = call(bridge, 'grep', [pattern: 'topsecret'])
        then: 'it IS searched - i.e. the sandbox guard, not the walk, is what excluded it above'
        allowed.count == 2
    }

    // =========================================================================
    // find
    // =========================================================================

    def 'find should match file names recursively'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'x'.bytes)
        Files.write(workDir.resolve('b.json'), 'x'.bytes)
        Files.createDirectories(workDir.resolve('sub'))
        Files.write(workDir.resolve('sub/c.txt'), 'x'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'find', [pattern: '*.txt'])
        then:
        result.count == 2
        result.paths.every { it.endsWith('.txt') }
        result.paths.any { it.endsWith('sub/c.txt') }
        result.truncated == false
    }

    def 'find should match a path-shaped glob against the relative path'() {
        given:
        def bridge = bridge()
        Files.createDirectories(workDir.resolve('results'))
        Files.write(workDir.resolve('results/stats.json'), 'x'.bytes)
        Files.write(workDir.resolve('stats.json'), 'x'.bytes)
        sandbox()
        when:
        def result = call(bridge, 'find', [pattern: 'results/*.json'])
        then:
        result.count == 1
        result.paths[0].endsWith('results/stats.json')
    }

    def 'find should filter by type'() {
        given:
        def bridge = bridge()
        Files.createDirectories(workDir.resolve('data'))
        Files.write(workDir.resolve('data.txt'), 'x'.bytes)
        sandbox()
        expect:
        call(bridge, 'find', [pattern: 'data*', type: 'dir']).count == 1
        call(bridge, 'find', [pattern: 'data*', type: 'file']).count == 1
        call(bridge, 'find', [pattern: 'data*']).count == 2
        call(bridge, 'find', [pattern: 'data*', type: 'socket']).error.contains('unknown type')
    }

    def 'find should cap the results and say so'() {
        given:
        def bridge = bridge()
        (1..10).each { Files.write(workDir.resolve("f${it}.txt"), 'x'.bytes) }
        sandbox()
        when:
        def result = call(bridge, 'find', [pattern: '*.txt', max_results: 4])
        then:
        result.count == 4
        result.truncated == true
        result.limit == 4
        result.truncated_reason == FilesystemTools.TRUNCATED_RESULTS
        when:
        def full = call(bridge, 'find', [pattern: '*.txt'])
        then:
        full.count == 10
        full.truncated == false
        full.truncated_reason == null
    }

    def 'find should cap the walk depth'() {
        given:
        def bridge = bridge()
        Files.createDirectories(workDir.resolve('l1/l2'))
        Files.write(workDir.resolve('top.txt'), 'x'.bytes)
        Files.write(workDir.resolve('l1/l2/deep.txt'), 'x'.bytes)
        sandbox()
        expect:
        call(bridge, 'find', [pattern: '*.txt', max_depth: 1]).count == 1
        call(bridge, 'find', [pattern: '*.txt']).count == 2
    }

    def 'find should refuse a root outside the sandbox and skip an escaping symlink'() {
        given:
        def bridge = bridge()
        escapingLink('secret.txt', outsideFile('secret.txt', 'topsecret'))
        Files.write(workDir.resolve('kept.txt'), 'x'.bytes)
        sandbox()
        when:
        def escaped = call(bridge, 'find', [pattern: '*.txt', path: '..'])
        then:
        escaped.error.contains('outside sandbox')
        when:
        def result = call(bridge, 'find', [pattern: '*.txt'])
        then: 'the symlink escaping the sandbox is not reported'
        result.paths.size() == 1
        result.paths[0].endsWith('kept.txt')
    }

    def 'find should reject a missing dir, a file and a bad glob'() {
        given:
        def bridge = bridge()
        Files.write(workDir.resolve('a.txt'), 'x'.bytes)
        sandbox()
        expect:
        call(bridge, 'find', [pattern: '*.txt', path: 'nope']).error.contains('path not found')
        call(bridge, 'find', [pattern: '*.txt', path: 'a.txt']).error.contains('not a directory')
        call(bridge, 'find', [pattern: '[bad']).error.contains('invalid glob')
        call(bridge, 'find', [:]).error.contains('missing required argument: pattern')
    }
}
