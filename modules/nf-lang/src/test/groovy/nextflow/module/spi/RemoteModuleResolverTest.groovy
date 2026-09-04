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

package nextflow.module.spi

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for the base directory rules of remote module include resolution.
 *
 * @author Ben Sherman <ben.sherman@seqera.io>
 */
class RemoteModuleResolverTest extends Specification {

    @TempDir
    Path projectDir

    def 'should resolve an include in a module script against the module directory'() {
        given:
        def moduleDir = projectDir.resolve('modules/nf-core/my-workflow')
        def nestedDir = moduleDir.resolve('modules/nf-core/star/align')

        expect:
        RemoteModuleResolver.resolveBaseDir(moduleDir.resolve('main.nf'), projectDir) == moduleDir
        and: 'at any depth'
        RemoteModuleResolver.resolveBaseDir(nestedDir.resolve('main.nf'), projectDir) == nestedDir
    }

    def 'should resolve an include in a non-module script against the project directory'() {
        given: 'an entry script in a subdirectory, as with a nested manifest.mainScript'
        def scriptDir = projectDir.resolve('workflows')

        expect:
        RemoteModuleResolver.resolveBaseDir(scriptDir.resolve('main.nf'), projectDir) == projectDir
        and: 'likewise for a script in the project root'
        RemoteModuleResolver.resolveBaseDir(projectDir.resolve('main.nf'), projectDir) == projectDir
    }

    def 'should fall back to the project directory when the including script is unknown'() {
        expect:
        RemoteModuleResolver.resolveBaseDir((Path) null, projectDir) == projectDir
        and:
        RemoteModuleResolver.resolveBaseDir((URI) null, projectDir) == projectDir
        and: 'a script outside the project is not treated as a vendored module'
        RemoteModuleResolver.resolveBaseDir(Path.of('/elsewhere/modules/a/b/main.nf'), projectDir) == projectDir
    }

    def 'should use the including script directory when the project directory is unknown'() {
        given: 'as when generating a spec for a module that has no spec yet'
        def moduleDir = projectDir.resolve('my-workflow')

        expect:
        RemoteModuleResolver.resolveBaseDir(moduleDir.resolve('main.nf'), null) == moduleDir
    }

    def 'should resolve a script identified by uri'() {
        given:
        def moduleDir = projectDir.resolve('modules/nf-core/my-workflow')
        Files.createDirectories(moduleDir)

        expect:
        RemoteModuleResolver.resolveBaseDir(moduleDir.resolve('main.nf').toUri(), projectDir) == moduleDir
        and: 'a uri that does not denote a file has no directory of its own'
        RemoteModuleResolver.resolveBaseDir(URI.create('untitled:Untitled-1'), projectDir) == projectDir
    }

}
