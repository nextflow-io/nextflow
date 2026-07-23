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

package nextflow.module

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for workflow-module validation in {@link ModuleValidator}.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleValidatorTest extends Specification {

    @TempDir
    Path tempDir

    private static final String PERMISSIVE_SCHEMA = '''\
        {
          "$schema": "https://json-schema.org/draft/2020-12/schema",
          "type": "object",
          "properties": { "name": {"type": "string"}, "description": {"type": "string"} },
          "required": ["name", "description"]
        }
        '''.stripIndent()

    private String schema() {
        final p = tempDir.resolve('schema.json')
        Files.writeString(p, PERMISSIVE_SCHEMA)
        return p.toString()
    }

    private Path moduleDir(String mainNf, String metaYml) {
        final d = Files.createDirectories(tempDir.resolve('mod'))
        Files.writeString(d.resolve('main.nf'), mainNf)
        Files.writeString(d.resolve('meta.yml'), metaYml)
        Files.writeString(d.resolve('README.md'), '# test module\n')
        return d
    }

    /**
     * Vendor a dependency under {@code <moduleDir>/modules/<scope>/<name>} at the given version,
     * mirroring nested per-module vendoring so {@code requires.modules} validation can find it.
     */
    private void vendorDependency(Path moduleDir, String scope, String name, String version) {
        final dep = moduleDir.resolve('modules').resolve(scope).resolve(name)
        Files.createDirectories(dep)
        Files.writeString(dep.resolve('main.nf'), "workflow DEP {\n take:\n x\n emit:\n x\n}\n")
        Files.writeString(dep.resolve('meta.yml'), "name: ${scope}/${name}\nversion: ${version}\nkind: Workflow\ndescription: vendored dep\n")
        Files.writeString(dep.resolve('README.md'), '# dep\n')
        ModuleInfo.save(dep, [checksum: ModuleChecksum.compute(dep), registryUrl: 'http://registry.com'])
    }

    def 'a workflow module that defines a workflow passes validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_in
                    main:
                    ch_out = ch_in
                    emit:
                    ch_out
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.isEmpty()
    }

    def 'a workflow module without a workflow definition fails validation' () {
        given:
        def dir = moduleDir(
            '''\
                process FOO {
                    """
                    echo hello
                    """
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.any { it.contains('must define a workflow') }
    }

    def 'a workflow module defining more than one workflow fails validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_in
                    emit:
                    ch_in
                }
                workflow BAR {
                    take:
                    ch_in
                    emit:
                    ch_in
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.any { it.contains('exactly one workflow') }
    }

    def 'a workflow meta.yml with matching input/output counts passes validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_a
                    ch_b
                    main:
                    ch_out = ch_a.mix(ch_b)
                    emit:
                    ch_out
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                input:
                  - name: ch_a
                    type: channel
                    description: first input
                  - name: ch_b
                    type: channel
                    description: second input
                output:
                  - name: ch_out
                    type: channel
                    description: the output
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.isEmpty()
    }

    def 'a workflow meta.yml with mismatched interface counts fails validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_a
                    ch_b
                    emit:
                    ch_a
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                input:
                  - name: ch_a
                    type: channel
                    description: only one declared, but the workflow takes two
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.any { it.contains('1 inputs but workflow declares 2 take') }
    }

    def 'a module whose requires.modules dependency is vendored at the pinned version passes' () {
        given:
        def dir = moduleDir(
            '''\
                include { DEP } from 'nf-core/dep'
                workflow FOO {
                    take:
                    ch_in
                    main:
                    ch_out = DEP(ch_in)
                    emit:
                    ch_out
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                requires:
                  modules:
                    - nf-core/dep@1.2.0
                '''.stripIndent())
        vendorDependency(dir, 'nf-core', 'dep', '1.2.0')

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.isEmpty()
    }

    def 'a module whose requires.modules dependency is not vendored fails validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_in
                    emit:
                    ch_in
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                requires:
                  modules:
                    - nf-core/dep@1.2.0
                '''.stripIndent())

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.any { it.contains('nf-core/dep@1.2.0') && it.contains('not vendored') }
    }

    def 'a module whose requires.modules dependency is vendored at a different version fails validation' () {
        given:
        def dir = moduleDir(
            '''\
                workflow FOO {
                    take:
                    ch_in
                    emit:
                    ch_in
                }
                '''.stripIndent(),
            '''\
                name: nf-core/demo_wf
                version: 1.0.0
                kind: Workflow
                description: a demo workflow module
                requires:
                  modules:
                    - nf-core/dep@2.0.0
                '''.stripIndent())
        vendorDependency(dir, 'nf-core', 'dep', '1.2.0')

        when:
        def errors = ModuleValidator.validate(dir, schema())

        then:
        errors.any { it.contains('vendored at version 1.2.0 but 2.0.0 is required') }
    }

}
