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
 * Unit tests for the {@code fs:} agent tools in {@link ModuleToolBridge}.
 *
 * <p>The tools are dispatched by their own wire names — {@code call('read', argsJson)},
 * {@code call('write', …)}, {@code call('ls', …)} — with no {@code operation} discriminator and
 * no aggregate {@code filesystem} tool: the family is six separately-named tools the model
 * selects between, so the tool name IS the operation. The bridge here is built with NO modules
 * and the whole {@code fs:} family selected.
 *
 * <p>The old {@code exists} operation has no successor and its cases are gone: {@code ls} and
 * {@code read} answer the same question, and a fourth way to ask it only cost the model a turn.
 *
 * <p>Scope note: the per-tool behaviour (argument handling, caps, the sandbox refusals for each
 * of the six) lives in {@link FilesystemToolsTest}. What is unique here is the interaction with
 * the MODULE side of the bridge — {@code whitelistOutputDirs} and its {@code isErrorResult}
 * guard, which is what lets a module output be read back without widening the sandbox.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleToolBridgeFilesystemTest extends Specification {

    @TempDir Path workDir

    def cleanup() {
        ModuleToolBridge.clearContext()
    }

    private ModuleToolBridge fsOnlyBridge() {
        // Build a bridge with NO modules and the whole `fs:` family selected
        new ModuleToolBridge(
            Collections.<WiredModuleTool>emptyList(),
            FilesystemTools.NAMES  // every fs: leaf
        )
    }

    private static Map parseJson(String json) {
        return (Map) new JsonSlurper().parseText(json)
    }

    private static String argsJson(Map args) {
        groovy.json.JsonOutput.toJson(args)
    }

    // -------------------------------------------------------------------------
    // filesystemDescriptors() carries the six fs tools when the family is selected,
    // and descriptors() -- the BROKERED half, i.e. `toolSpecs` -- never does (§5)
    // -------------------------------------------------------------------------

    def 'filesystem descriptors should carry the six fs tools under their own names when selected'() {
        given:
        def bridge = fsOnlyBridge()
        expect:
        (bridge.filesystemDescriptors()*.name as Set) == (['read','write','edit','ls','grep','find'] as Set)
        and: 'never an aggregate tool the model would have to discriminate with an argument'
        !bridge.filesystemDescriptors().any { it.name == 'filesystem' }
        and: 'and they stay OUT of the brokered descriptors that become `toolSpecs`'
        bridge.descriptors().isEmpty()
    }

    def 'descriptors should carry no fs tool when none is selected'() {
        given:
        def bridge = new ModuleToolBridge(
            Collections.<WiredModuleTool>emptyList(),
            Collections.<String>emptyList()
        )
        expect:
        bridge.descriptors().isEmpty()
        bridge.filesystemDescriptors().isEmpty()
    }

    // -------------------------------------------------------------------------
    // no context → error
    // -------------------------------------------------------------------------

    def 'an fs call without context returns error'() {
        given:
        def bridge = fsOnlyBridge()
        // no context set
        when:
        def result = parseJson(bridge.call('read', argsJson([path: 'test.txt'])))
        then:
        result.error != null
        result.error.contains('no sandbox context')
    }

    // -------------------------------------------------------------------------
    // write
    // -------------------------------------------------------------------------

    def 'write creates a file in workDir'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('write', argsJson([
            path: 'output.txt',
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
        def result = bridge.call('read', argsJson([path: 'data.txt']))
        then:
        // readOrHandle returns inline string for .txt files; the JSON wraps it under 'content'
        def parsed = parseJson(result)
        parsed.content == 'the content'
    }

    // -------------------------------------------------------------------------
    // ls
    // -------------------------------------------------------------------------

    def 'ls returns directory entries'() {
        given:
        def bridge = fsOnlyBridge()
        Files.write(workDir.resolve('a.txt'), 'a'.bytes)
        Files.write(workDir.resolve('b.txt'), 'b'.bytes)
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('ls', argsJson([path: '.'])))
        then:
        result.entries instanceof List
        (result.entries*.name as Set).containsAll(['a.txt', 'b.txt'])
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
        def result = parseJson(bridge.call('read', argsJson([path: '../secret.txt'])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    def 'write outside sandbox returns error'() {
        given:
        def bridge = fsOnlyBridge()
        ModuleToolBridge.setContext(new DispatchContext(workDir))
        when:
        def result = parseJson(bridge.call('write', argsJson([path: '../evil.txt', content: 'bad'])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    // -------------------------------------------------------------------------
    // addReadablePath allows reads from module-output paths OUTSIDE workDir
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
        ctx.addReadablePath(outsideDir)
        ModuleToolBridge.setContext(ctx)
        when:
        // (a) read of a file inside the whitelisted outside dir SUCCEEDS
        def result = parseJson(bridge.call('read', argsJson([
            path: outsideDir.resolve('result.txt').toAbsolutePath().toString()
        ])))
        then:
        result.content == 'module result'
    }

    def 'read from a non-whitelisted outside dir returns error'() {
        given:
        def bridge = fsOnlyBridge()
        // sandboxDir is the real workDir; outsideDir is a sibling NOT added to readablePaths
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('sandbox2-' + System.nanoTime()))
        def notWhitelisted = Files.createDirectories(tmp.resolve('not-whitelisted-' + System.nanoTime()))
        Files.write(notWhitelisted.resolve('secret.txt'), 'secret'.bytes)
        def ctx = new DispatchContext(sandboxDir)
        // do NOT add notWhitelisted to ctx
        ModuleToolBridge.setContext(ctx)
        when:
        // (b) read of a file in a DIFFERENT outside dir that is NOT whitelisted returns {"error":...}
        def result = parseJson(bridge.call('read', argsJson([
            path: notWhitelisted.resolve('secret.txt').toAbsolutePath().toString()
        ])))
        then:
        result.error != null
        result.error.contains('outside sandbox')
    }

    // -------------------------------------------------------------------------
    // I2: whitelistOutputDirs auto-adds the produced file so an fs read succeeds
    // -------------------------------------------------------------------------

    /**
     * Exercises the auto-whitelist path: calling {@link ModuleToolBridge#whitelistOutputDirs}
     * with a parsed result Map that contains an absolute file path OUTSIDE the sandbox work dir
     * must add THAT FILE — not its parent — to the context's readablePaths, so a subsequent
     * {@code read} of that path succeeds while a SIBLING of the output stays rejected.
     * Whitelisting the parent would grant every sibling too, and for an output landing outside
     * the work tree that is a directory of content no cache key covers.
     *
     * Approach chosen: focused unit test calling the package-visible {@code whitelistOutputDirs}
     * directly (rather than wiring a live process). This keeps the test lightweight and precisely
     * exercises the auto-whitelisting logic without requiring a full Nextflow process harness.
     * The {@code read} call through {@link ModuleToolBridge#call} then verifies that the sandbox
     * guard correctly allows the whitelisted path and rejects the non-whitelisted one.
     */
    def 'whitelistOutputDirs adds the produced file only, not its parent dir'() {
        given: 'a sandbox work dir and a module output dir OUTSIDE the sandbox'
        def bridge = fsOnlyBridge()
        def tmp = workDir.getParent()
        def sandboxDir = Files.createDirectories(tmp.resolve('wl-sandbox-' + System.nanoTime()))
        def moduleOutDir = Files.createDirectories(tmp.resolve('wl-module-out-' + System.nanoTime()))
        def nonWhitelistedDir = Files.createDirectories(tmp.resolve('wl-other-' + System.nanoTime()))

        // write the module output file outside the sandbox
        def outputFile = moduleOutDir.resolve('result.txt')
        Files.write(outputFile, 'module output content'.bytes)
        // a SIBLING of the output, in the very same dir -- the file the old parent-dir
        // whitelist would have exposed
        def siblingFile = moduleOutDir.resolve('sibling-secret.txt')
        Files.write(siblingFile, 'not an output'.bytes)
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

        then: 'the produced FILE is whitelisted, and its parent dir is not'
        ctx.readablePaths.contains(outputFile)
        !ctx.readablePaths.contains(moduleOutDir)

        and: 'read of the output file SUCCEEDS (auto-whitelisted)'
        def readResult = parseJson(bridge.call('read', argsJson([
            path: outputFile.toAbsolutePath().toString()
        ])))
        readResult.content == 'module output content'

        and: 'a SIBLING of the output, in the same dir, is REJECTED'
        def siblingResult = parseJson(bridge.call('read', argsJson([
            path: siblingFile.toAbsolutePath().toString()
        ])))
        siblingResult.error != null
        siblingResult.error.contains('outside sandbox')

        and: 'read of a file in a non-whitelisted sibling dir is REJECTED'
        def rejectResult = parseJson(bridge.call('read', argsJson([
            path: nonWhitelistedDir.resolve('secret.txt').toAbsolutePath().toString()
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

    def 'whitelistOutputDirs auto-adds the produced path of non-error results'() {
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

        then: 'the produced file is now in readablePaths'
        ctx.readablePaths.contains(outputFile)

        and: 'read of the output file SUCCEEDS (auto-whitelisted)'
        def readResult = parseJson(bridge.call('read', argsJson([
            path: outputFile.toAbsolutePath().toString()
        ])))
        readResult.content == 'module output content'
    }
}
