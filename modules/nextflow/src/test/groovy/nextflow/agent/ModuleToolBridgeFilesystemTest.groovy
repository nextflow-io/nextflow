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
    // addReadableDir allows reads from module-output paths
    // -------------------------------------------------------------------------

    def 'read from a whitelisted readable dir succeeds'() {
        given:
        def bridge = fsOnlyBridge()
        def moduleOutputDir = Files.createDirectories(workDir.resolve('module-out'))
        Files.write(moduleOutputDir.resolve('result.txt'), 'module result'.bytes)
        def ctx = new DispatchContext(workDir)
        ctx.addReadableDir(moduleOutputDir)
        ModuleToolBridge.setContext(ctx)
        when:
        // Use absolute path for a file in module-output dir (not inside workDir itself)
        def result = parseJson(bridge.call('filesystem', argsJson([
            path: moduleOutputDir.resolve('result.txt').toAbsolutePath().toString(),
            operation: 'read'
        ])))
        then:
        result.content == 'module result'
    }
}
