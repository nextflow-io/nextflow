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
 * Unit tests for the {@code filesystem} agent tool path in {@link ModuleToolBridge}.
 * These tests exercise the {@code call('filesystem', argsJson)} dispatch without
 * any real dataflow processes - the bridge is constructed with no modules and only
 * the {@code filesystem} capability flag enabled.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleToolBridgeFilesystemTest extends Specification {

    @TempDir Path workDir

    def cleanup() {
        ModuleToolBridge.clearContext()
    }

    private ModuleToolBridge fsOnlyBridge() {
        // Build a bridge with NO modules but filesystem capability enabled
        new ModuleToolBridge(
            Collections.<String, nextflow.script.ProcessDef>emptyMap(),
            Collections.<String, nextflow.module.ModuleSpec>emptyMap(),
            Collections.<String, ModuleToolBridge.RegistryMeta>emptyMap(),
            true  // filesystemEnabled
        )
    }

    private static Map parseJson(String json) {
        return (Map) new JsonSlurper().parseText(json)
    }

    private static String argsJson(Map args) {
        groovy.json.JsonOutput.toJson(args)
    }

    // -------------------------------------------------------------------------
    // descriptor() includes the filesystem tool when enabled
    // -------------------------------------------------------------------------

    def 'descriptors should include filesystem tool when filesystemEnabled'() {
        given:
        def bridge = fsOnlyBridge()
        expect:
        bridge.descriptors().any { it.name == 'filesystem' }
    }

    def 'descriptors should NOT include filesystem tool when not enabled'() {
        given:
        def bridge = new ModuleToolBridge(
            Collections.<String, nextflow.script.ProcessDef>emptyMap(),
            Collections.<String, nextflow.module.ModuleSpec>emptyMap(),
            Collections.<String, ModuleToolBridge.RegistryMeta>emptyMap(),
            false
        )
        expect:
        !bridge.descriptors().any { it.name == 'filesystem' }
    }

    // -------------------------------------------------------------------------
    // no context → error
    // -------------------------------------------------------------------------

    def 'filesystem call without context returns error'() {
        given:
        def bridge = fsOnlyBridge()
        // no context set
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([path: 'test.txt', operation: 'exists'])))
        then:
        result.error != null
        result.error.contains('no sandbox context')
    }

    // -------------------------------------------------------------------------
    // exists
    // -------------------------------------------------------------------------

    def 'exists returns true for an existing file in workDir'() {
        given:
        def bridge = fsOnlyBridge()
        def file = Files.write(workDir.resolve('hello.txt'), 'world'.bytes)
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([path: 'hello.txt', operation: 'exists'])))
        then:
        result.exists == true
    }

    def 'exists returns false for a non-existing file in workDir'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([path: 'missing.txt', operation: 'exists'])))
        then:
        result.exists == false
    }

    // -------------------------------------------------------------------------
    // write
    // -------------------------------------------------------------------------

    def 'write creates a file in workDir'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([
            path: 'output.txt',
            operation: 'write',
            content: 'hello world'
        ])))
        then:
        result.error == null
        result.bytes == 11
        workDir.resolve('output.txt').text == 'hello world'
    }

    // -------------------------------------------------------------------------
    // read
    // -------------------------------------------------------------------------

    def 'read returns the content of a file in workDir'() {
        given:
        def bridge = fsOnlyBridge()
        Files.write(workDir.resolve('data.txt'), 'the content'.bytes)
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = bridge.call('filesystem', argsJson([path: 'data.txt', operation: 'read']))
        then:
        // readOrHandle returns inline string for .txt files; the JSON wraps it under 'content'
        def parsed = parseJson(result)
        parsed.content == 'the content'
    }

    // -------------------------------------------------------------------------
    // list
    // -------------------------------------------------------------------------

    def 'list returns directory entries'() {
        given:
        def bridge = fsOnlyBridge()
        Files.write(workDir.resolve('a.txt'), 'a'.bytes)
        Files.write(workDir.resolve('b.txt'), 'b'.bytes)
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([path: '.', operation: 'list'])))
        then:
        result.entries instanceof List
        result.entries.containsAll(['a.txt', 'b.txt'])
    }

    // -------------------------------------------------------------------------
    // sandbox enforcement
    // -------------------------------------------------------------------------

    def 'read of a path outside the sandbox returns error'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        // Use a path that escapes the workDir via ..
        def result = parseJson(bridge.call('filesystem', argsJson([path: '../secret.txt', operation: 'read'])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    def 'write outside sandbox returns error'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('filesystem', argsJson([path: '../evil.txt', operation: 'write', content: 'bad'])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    // -------------------------------------------------------------------------
    // addReadableDir allows reads from module-output paths OUTSIDE workDir
    // -------------------------------------------------------------------------

    def 'read from a whitelisted readable dir outside workDir succeeds'() {
        given:
        def bridge = fsOnlyBridge()
        // Use a SIBLING of workDir so module-out-sibling is NOT inside workDir;
        // this exercises the whitelist branch (not the isInside(workDir) branch).
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('sandbox-' + System.nanoTime()))
        def outsideDir = Files.createDirectories(tmp.resolve('module-out-sibling-' + System.nanoTime()))
        Files.write(outsideDir.resolve('result.txt'), 'module result'.bytes)
        def ctx = new DispatchContext(sandboxDir)
        ctx.addReadableDir(outsideDir)
        ModuleToolBridge.setContext(ctx)
        when:
        // (a) read of a file inside the whitelisted outside dir SUCCEEDS
        def result = parseJson(bridge.call('filesystem', argsJson([
            path: outsideDir.resolve('result.txt').toAbsolutePath().toString(),
            operation: 'read'
        ])))
        then:
        result.content == 'module result'
    }

    def 'read from a non-whitelisted outside dir returns error'() {
        given:
        def bridge = fsOnlyBridge()
        // sandboxDir is the real workDir; outsideDir is a sibling NOT added to readableDirs
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('sandbox2-' + System.nanoTime()))
        def notWhitelisted = Files.createDirectories(tmp.resolve('not-whitelisted-' + System.nanoTime()))
        Files.write(notWhitelisted.resolve('secret.txt'), 'secret'.bytes)
        def ctx = new DispatchContext(sandboxDir)
        // do NOT add notWhitelisted to ctx
        ModuleToolBridge.setContext(ctx)
        when:
        // (b) read of a file in a DIFFERENT outside dir that is NOT whitelisted returns {"error":...}
        def result = parseJson(bridge.call('filesystem', argsJson([
            path: notWhitelisted.resolve('secret.txt').toAbsolutePath().toString(),
            operation: 'read'
        ])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    // -------------------------------------------------------------------------
    // I2: whitelistOutputDirs auto-adds parent dir so filesystem read succeeds
    // -------------------------------------------------------------------------

    /**
     * Exercises the auto-whitelist path: calling {@link ModuleToolBridge#whitelistOutputDirs}
     * with a parsed result Map that contains an absolute file path OUTSIDE the sandbox work dir
     * must add that file's parent dir to the context's readableDirs, so a subsequent
     * {@code filesystem} {@code read} of that path succeeds — whereas a sibling path in a
     * non-whitelisted dir is still rejected.
     *
     * Approach chosen: focused unit test calling the package-visible {@code whitelistOutputDirs}
     * directly (rather than wiring a live process). This keeps the test lightweight and precisely
     * exercises the auto-whitelisting logic without requiring a full Nextflow process harness.
     * The filesystem {@code read} call through {@link ModuleToolBridge#call} then verifies that
     * the sandbox guard correctly allows the whitelisted path and rejects the non-whitelisted one.
     */
    def 'whitelistOutputDirs auto-adds parent dir enabling filesystem read of module output'() {
        given: 'a sandbox work dir and a module output dir OUTSIDE the sandbox'
        def bridge = fsOnlyBridge()
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('wl-sandbox-' + System.nanoTime()))
        def moduleOutDir = Files.createDirectories(tmp.resolve('wl-module-out-' + System.nanoTime()))
        def nonWhitelistedDir = Files.createDirectories(tmp.resolve('wl-other-' + System.nanoTime()))

        // write the module output file outside the sandbox
        def outputFile = moduleOutDir.resolve('result.txt')
        Files.write(outputFile, 'module output content'.bytes)
        // write a file in the non-whitelisted sibling dir
        Files.write(nonWhitelistedDir.resolve('secret.txt'), 'secret'.bytes)

        // set up a dispatch context with the sandbox dir
        def ctx = new DispatchContext(sandboxDir)
        ModuleToolBridge.setContext(ctx)

        when: 'whitelistOutputDirs is called with a parsed result map containing the output file path'
        // simulate what call() does after a successful module-tool dispatch: the parsed result
        // contains the absolute path of the output file
        def resultMap = [output: outputFile.toAbsolutePath().toString()]
        ModuleToolBridge.whitelistOutputDirs(resultMap)

        then: 'the output file parent dir is now in readableDirs'
        ctx.readableDirs.contains(moduleOutDir)

        and: 'filesystem read of the output file SUCCEEDS (auto-whitelisted)'
        def readResult = parseJson(bridge.call('filesystem', argsJson([
            path: outputFile.toAbsolutePath().toString(),
            operation: 'read'
        ])))
        readResult.content == 'module output content'

        and: 'filesystem read of a file in a non-whitelisted sibling dir is REJECTED'
        def rejectResult = parseJson(bridge.call('filesystem', argsJson([
            path: nonWhitelistedDir.resolve('secret.txt').toAbsolutePath().toString(),
            operation: 'read'
        ])))
        rejectResult.error != null
        rejectResult.error.contains('outside sandbox')
    }

    // -------------------------------------------------------------------------
    // I1: guard predicate isErrorResult detects error vs. non-error results
    // -------------------------------------------------------------------------

    /**
     * Tests the extracted guard predicate {@link ModuleToolBridge#isErrorResult}.
     * This predicate is used in {@link ModuleToolBridge#call} to skip
     * {@code whitelistOutputDirs} when the result is an error, preventing absolute
     * paths in error messages from widening the filesystem sandbox.
     */
    def 'isErrorResult detects error-shaped results (Map with error key)'() {
        given: 'parse an error result JSON containing an absolute path in the error message'
        def errorJson = '{"error":"module failed at /etc/secret/data.txt"}'
        def parsed = new JsonSlurper().parseText(errorJson)

        expect: 'isErrorResult returns true for error-shaped results'
        ModuleToolBridge.isErrorResult(parsed) == true

        and: 'whitelistOutputDirs must be skipped for error results (guarded by isErrorResult)'
        // this is what call() checks; the guard prevents whitelistOutputDirs
        // from being called on error results, so paths in error messages stay out of the whitelist
    }

    def 'isErrorResult returns false for non-error results'() {
        given: 'parse a normal (non-error) result JSON with a file path'
        def normalJson = '{"out":"/tmp/work/abc/result.fa"}'
        def parsed = new JsonSlurper().parseText(normalJson)

        expect: 'isErrorResult returns false for non-error results'
        ModuleToolBridge.isErrorResult(parsed) == false

        and: 'whitelistOutputDirs proceeds for non-error results (guarded by !isErrorResult)'
        // this is what call() checks; when the result is NOT an error,
        // whitelistOutputDirs is called to add output file paths to the whitelist
    }

    def 'whitelistOutputDirs auto-adds parent dir of non-error paths'() {
        given: 'a sandbox work dir and a module output dir OUTSIDE the sandbox'
        def bridge = fsOnlyBridge()
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('wl-sandbox-' + System.nanoTime()))
        def moduleOutDir = Files.createDirectories(tmp.resolve('wl-module-out-' + System.nanoTime()))

        // write the module output file outside the sandbox
        def outputFile = moduleOutDir.resolve('result.txt')
        Files.write(outputFile, 'module output content'.bytes)

        // set up a dispatch context with the sandbox dir
        def ctx = new DispatchContext(sandboxDir)
        ModuleToolBridge.setContext(ctx)

        when: 'a non-error result JSON containing the output file path is parsed and whitelisted'
        def resultJson = JsonOutput.toJson([output: outputFile.toAbsolutePath().toString()])
        def resultParsed = new JsonSlurper().parseText(resultJson)
        // guard check: the result is NOT an error, so whitelistOutputDirs IS called
        if( !ModuleToolBridge.isErrorResult(resultParsed) )
            ModuleToolBridge.whitelistOutputDirs(resultParsed)

        then: 'the output file parent dir is now in readableDirs'
        ctx.readableDirs.contains(moduleOutDir)

        and: 'filesystem read of the output file SUCCEEDS (auto-whitelisted)'
        def readResult = parseJson(bridge.call('filesystem', argsJson([
            path: outputFile.toAbsolutePath().toString(),
            operation: 'read'
        ])))
        readResult.content == 'module output content'
    }
}
