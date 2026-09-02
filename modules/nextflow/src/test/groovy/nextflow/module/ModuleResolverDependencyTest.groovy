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

import io.seqera.npr.api.schema.v1.ListModuleReleasesResponse
import io.seqera.npr.api.schema.v1.Module
import io.seqera.npr.api.schema.v1.ModuleRelease
import io.seqera.npr.client.RegistryClient
import nextflow.exception.AbortOperationException
import nextflow.util.VersionNumber
import org.yaml.snakeyaml.Yaml
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for nested transitive dependency vendoring in
 * {@link ModuleResolver#installWithDependencies} (ADR v2, PR #7342).
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleResolverDependencyTest extends Specification {

    @TempDir
    Path tempDir

    // module universe: fullName -> ( version -> requires.modules list )
    private Map<String, Map<String, List<String>>> universe = [:]

    private void module(String name, String version, List<String> requires = []) {
        universe.computeIfAbsent(name, { [:] }).put(version, requires)
    }

    private void buildBundle(String name, String version, List<String> requires, Path dest) {
        final dir = Files.createTempDirectory('mod')
        dir.resolve('main.nf').text = "workflow FOO {\n}\n"
        final meta = [name: name, version: version, kind: 'Workflow', description: "test module ${name}".toString(),
                      requires: [nextflow: '>=24.04.0']]
        if( requires )
            meta.requires.modules = requires
        dir.resolve('meta.yml').text = new Yaml().dump(meta)
        ModuleStorage.createBundle(dir, dest)
    }

    private RegistryClient mockClient() {
        final client = Mock(RegistryClient)
        client.listModuleReleases(_) >> { String n ->
            final resp = new ListModuleReleasesResponse()
            resp.releases = universe.getOrDefault(n, [:]).keySet().collect { v -> new ModuleRelease().version(v) }
            return resp
        }
        client.getModule(_) >> { String n ->
            final versions = new ArrayList<String>(universe.getOrDefault(n, [:]).keySet())
            final latest = versions.max { a, b -> new VersionNumber(a) <=> new VersionNumber(b) }
            return new Module().latest(new ModuleRelease().version(latest))
        }
        client.downloadModuleRelease(_, _, _) >> { String n, String v, Path dest ->
            buildBundle(n, v, universe[n][v], dest)
            return "oci://${n}:${v}"
        }
        return client
    }

    private boolean installedAt(String relDir) {
        return Files.exists(tempDir.resolve(relDir).resolve('main.nf'))
    }

    private String versionAt(String relDir) {
        return ModuleSpecFactory.fromYaml(tempDir.resolve(relDir).resolve('meta.yml')).version
    }

    // simulate a user hand-editing an installed module's meta.yml requires.modules
    private void editRequires(String relDir, List<String> requires) {
        final metaFile = tempDir.resolve(relDir).resolve('meta.yml')
        final meta = new Yaml().load(metaFile.text) as Map
        final req = (meta.requires ?: [:]) as Map
        if( requires )
            req.modules = requires
        else
            req.remove('modules')
        meta.requires = req
        metaFile.text = new Yaml().dump(meta)
    }

    def 'a freshly installed workflow module with dependencies reads as VALID' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        def storage = new ModuleStorage(tempDir)

        when:
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')

        then: 'the checksum covers the vendored subtree, so no false MODIFIED status'
        storage.getInstalledModule(ModuleReference.parse('nf-core/aln')).integrity == ModuleIntegrity.VALID
        and: 'the vendored dependency is itself VALID'
        new ModuleStorage(tempDir.resolve('modules/nf-core/aln'))
            .getInstalledModule(ModuleReference.parse('nf-core/samtools/sort')).integrity == ModuleIntegrity.VALID
    }

    def 'editing an installed module own file marks it MODIFIED' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        def storage = new ModuleStorage(tempDir)
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')

        when:
        def mainNf = tempDir.resolve('modules/nf-core/aln/main.nf')
        mainNf.text = mainNf.text + "\n// local edit\n"

        then:
        storage.getInstalledModule(ModuleReference.parse('nf-core/aln')).integrity == ModuleIntegrity.MODIFIED
    }

    def 'editing a vendored dependency by hand marks the parent MODIFIED' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        def storage = new ModuleStorage(tempDir)
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')

        when: 'a vendored dependency file is edited by hand'
        def depMain = tempDir.resolve('modules/nf-core/aln/modules/nf-core/samtools/sort/main.nf')
        depMain.text = depMain.text + "\n// hand edit of a vendored dep\n"

        then: 'the change is surfaced through the parent module integrity (Option 1: subtree checksum)'
        storage.getInstalledModule(ModuleReference.parse('nf-core/aln')).integrity == ModuleIntegrity.MODIFIED
    }

    def 'update-deps installs a newly declared dependency' () {
        given:
        module('nf-core/aln', '1.0.0', [])            // parent initially has no dependencies
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')
        assert !installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')

        when: 'the user declares a dependency in meta.yml and updates deps'
        editRequires('modules/nf-core/aln', ['nf-core/samtools/sort@1.0.0'])
        resolver.updateDependencies(ModuleReference.parse('nf-core/aln'))

        then:
        installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')
    }

    def 'update-deps updates a dependency to the version declared in meta.yml' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        module('nf-core/samtools/sort', '2.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')
        assert versionAt('modules/nf-core/aln/modules/nf-core/samtools/sort') == '1.0.0'

        when:
        editRequires('modules/nf-core/aln', ['nf-core/samtools/sort@2.0.0'])
        resolver.updateDependencies(ModuleReference.parse('nf-core/aln'))

        then:
        versionAt('modules/nf-core/aln/modules/nf-core/samtools/sort') == '2.0.0'
    }

    def 'update-deps prunes a dependency removed from meta.yml' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')
        assert installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')

        when: 'the dependency is removed from meta.yml and deps are updated'
        editRequires('modules/nf-core/aln', [])
        resolver.updateDependencies(ModuleReference.parse('nf-core/aln'))

        then:
        !installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')
    }

    def 'update-deps errors and keeps a removed dependency that has local modifications' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')
        and: 'the vendored dependency is hand-edited'
        def depMain = tempDir.resolve('modules/nf-core/aln/modules/nf-core/samtools/sort/main.nf')
        depMain.text = depMain.text + "\n// local edit\n"

        when: 'the dependency is removed from meta.yml and deps are updated'
        editRequires('modules/nf-core/aln', [])
        resolver.updateDependencies(ModuleReference.parse('nf-core/aln'))

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('local modifications')
        and: 'the modified orphan is not removed'
        installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')
    }

    def 'update-deps leaves the parent module untouched (stays MODIFIED)' () {
        given:
        module('nf-core/aln', '1.0.0', [])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())
        def storage = new ModuleStorage(tempDir)
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')

        when: 'the user edits the parent meta.yml to add a dependency, then updates deps'
        editRequires('modules/nf-core/aln', ['nf-core/samtools/sort@1.0.0'])
        resolver.updateDependencies(ModuleReference.parse('nf-core/aln'))

        then: 'the dependency is vendored'
        installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')
        and: 'the parent stays MODIFIED (unpublished meta.yml edit; checksum not refreshed)'
        storage.getInstalledModule(ModuleReference.parse('nf-core/aln')).integrity == ModuleIntegrity.MODIFIED
    }

    def 'should install a module and vendor its dependencies under its own nested modules/ dir' () {
        given:
        module('nf-core/aln', '1.0.0', ['nf-core/samtools/sort@1.0.0'])
        module('nf-core/samtools/sort', '1.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())

        when:
        resolver.installWithDependencies(ModuleReference.parse('nf-core/aln'), '1.0.0')

        then:
        installedAt('modules/nf-core/aln')
        // dependency is vendored under the module's OWN nested modules/ directory
        installedAt('modules/nf-core/aln/modules/nf-core/samtools/sort')
    }

    def 'should vendor a shared dependency independently per consumer (duplication, no flattening)' () {
        given:
        module('nf-core/top', '1.0.0', ['nf-core/a@1.0.0', 'nf-core/b@1.0.0'])
        module('nf-core/a', '1.0.0', ['nf-core/c@1.0.0'])
        module('nf-core/b', '1.0.0', ['nf-core/c@2.0.0'])  // different version -- no conflict under nested model
        module('nf-core/c', '1.0.0')
        module('nf-core/c', '2.0.0')
        def resolver = new ModuleResolver(tempDir, mockClient())

        when:
        resolver.installWithDependencies(ModuleReference.parse('nf-core/top'), '1.0.0')

        then:
        // c is vendored twice -- once under a, once under b -- at each pinned version
        versionAt('modules/nf-core/top/modules/nf-core/a/modules/nf-core/c') == '1.0.0'
        versionAt('modules/nf-core/top/modules/nf-core/b/modules/nf-core/c') == '2.0.0'
    }

    def 'should detect a dependency cycle' () {
        given:
        module('nf-core/a', '1.0.0', ['nf-core/b@1.0.0'])
        module('nf-core/b', '1.0.0', ['nf-core/a@1.0.0'])
        def resolver = new ModuleResolver(tempDir, mockClient())

        when:
        resolver.installWithDependencies(ModuleReference.parse('nf-core/a'), '1.0.0')

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('cycle')
    }

    def 'resolve with autoInstall installs a workflow module and its nested deps' () {
        given:
        module('nf-core/mafft_align', '0.0.0-test001', ['nf-core/mafft/align@0.0.0-test001'])
        module('nf-core/mafft/align', '0.0.0-test001')
        def resolver = new ModuleResolver(tempDir, mockClient())

        when:
        // this is the include-resolution path (autoInstall = true, version = null)
        def main = resolver.resolve(ModuleReference.parse('nf-core/mafft_align'), null, true)

        then:
        installedAt('modules/nf-core/mafft_align')
        installedAt('modules/nf-core/mafft_align/modules/nf-core/mafft/align')
        main.toString().endsWith('modules/nf-core/mafft_align/main.nf')
    }

}
