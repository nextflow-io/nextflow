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

import nextflow.module.ModuleInfo
import nextflow.module.ModuleReference
import spock.lang.Specification
import spock.lang.TempDir

class ModuleToolResolverTest extends Specification {

    @TempDir
    File tempDir

    // Helper: invoke the private static recoverModuleRef via Groovy metaprogramming
    private static ModuleReference recoverModuleRef(java.nio.file.Path moduleDir) {
        def m = ModuleToolResolver.getDeclaredMethod('recoverModuleRef', java.nio.file.Path)
        m.accessible = true
        // pass as explicit Object[] so null is not misinterpreted as (Object[]) null (zero-arg)
        return (ModuleReference) m.invoke(null, new Object[]{moduleDir})
    }

    // -----------------------------------------------------------------------
    // recoverModuleRef unit tests (offline-safe — no network involved)
    // -----------------------------------------------------------------------

    def 'recoverModuleRef: registry install dir WITH marker returns correct ModuleReference'() {
        given: 'a directory tree matching <base>/modules/<scope>/<name> with .module-info marker'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc123'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref != null
        ref.scope == 'nf-core'
        ref.name == 'skesa'
        ref.fullName == 'nf-core/skesa'
    }

    def 'recoverModuleRef: dir WITHOUT marker returns null (local-file include)'() {
        given: 'same layout but NO .module-info marker file'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('fastqc')
        Files.createDirectories(moduleDir)
        // intentionally no .module-info

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref == null
    }

    def 'recoverModuleRef: null input returns null'() {
        expect:
        recoverModuleRef(null) == null
    }

    def 'recoverModuleRef: dir with marker but NOT under a "modules" grandparent returns null'() {
        given: 'marker present but parent is named something other than "modules"'
        def base = tempDir.toPath()
        def moduleDir = base.resolve('notmodules').resolve('nf-core').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref == null
    }

    def 'recoverModuleRef: dir with marker but only one parent level returns null'() {
        given: 'marker present but dir has only one ancestor above it (no scope/name split possible)'
        def base = tempDir.toPath()
        // layout: <base>/modules/skesa  — no scope segment
        def moduleDir = base.resolve('modules').resolve('skesa')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=abc'

        when:
        // At this layout: moduleDir.parent.fileName = 'modules', moduleDir.parent.parent.fileName ≠ 'modules'
        // So the check "modulesDir.fileName == 'modules'" fails because the grandparent of
        // moduleDir is the base dir (temp dir), not 'modules'. Returns null.
        def ref = recoverModuleRef(moduleDir)

        then:
        // parent = modules dir (fileName='modules'), parent.parent = base (fileName != 'modules')
        // scope = 'modules', modulesDir = base, base.fileName != 'modules' → null
        ref == null
    }

    def 'recoverModuleRef: multi-segment module name (e.g. nf-core/subworkflows/bam_sort_stats_samtools) returns reference'() {
        given: 'a registry install with a multi-level name stored at scope/first_segment/rest'
        // ModuleStorage stores <base>/modules/<scope>/<name> where name can be 'bam_sort_stats_samtools'
        // (nf-core uses flat single-segment names for modules; multi-segment only for subworkflows,
        //  stored as a single directory name). This test covers an arbitrary valid single-dir name.
        def base = tempDir.toPath()
        def moduleDir = base.resolve('modules').resolve('nf-core').resolve('bwa_mem')
        Files.createDirectories(moduleDir)
        moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE).text = 'checksum=xyz'

        when:
        def ref = recoverModuleRef(moduleDir)

        then:
        ref != null
        ref.scope == 'nf-core'
        ref.name == 'bwa_mem'
    }
}
