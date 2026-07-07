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

package nextflow.install2r

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification

class Install2rCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should create install2r env prefix path for a string env' () {
        given:
        def ENV = 'dplyr ggplot2'
        def cache = Spy(Install2rCache)
        def BASE = Paths.get('/i2r/envs')

        when:
        def prefix = cache.install2rPrefixPath(ENV)
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/i2r/envs/env-')
    }

    def 'should create the correct install2.r command for package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(Install2rCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalInstall2rEnv0('dplyr ggplot2', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('install2.r') && c.contains('-l') && c.contains('dplyr ggplot2')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should build a command with --error and both packages for a multi-package env' () {
        // second valid-config shape: a multi-package list must produce a single
        // install2.r invocation carrying --error, the -l library flag, and BOTH
        // package names. install2.r has no manifest-file concept, so packages are
        // always passed inline as bare names on the command line.
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-multi')
        def cache = Spy(Install2rCache)
        cache.@installOptions = null
        cache.@repos = 'https://cloud.r-project.org'
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalInstall2rEnv0('dplyr ggplot2', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('install2.r') &&
                c.contains('--error') &&
                c.contains('-l') &&
                c.contains('dplyr') &&
                c.contains('ggplot2')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should treat a slash-containing spec that is not a directory as a package name' () {
        // Edge case guarded by install2rPrefixPath: only an EXISTING directory is
        // returned verbatim. A path-with-slash that does not resolve to a directory
        // (here: an install2.r-style repo-qualified package reference) is hashed
        // into a normal cache env path instead of being used as a literal library.
        given:
        def cache = Spy(Install2rCache)
        def BASE = Paths.get('/i2r/envs')

        when:
        def prefix = cache.install2rPrefixPath('some/not/a/dir')
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/i2r/envs/env-')
    }

    def 'should return a real existing directory spec verbatim as the library path' () {
        // Complement to the edge case above: when the slash-containing spec IS an
        // existing directory it is used directly as the library prefix (no hashing,
        // no getCacheDir lookup).
        given:
        def folder = Files.createTempDirectory('lib')
        def cache = Spy(Install2rCache)

        when:
        def prefix = cache.install2rPrefixPath(folder.toString())
        then:
        0 * cache.getCacheDir()
        prefix.toString() == folder.toString()

        cleanup:
        folder?.deleteDir()
    }

    def 'should apply a per-process install-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(Install2rCache)
        cache.@installOptions = '--from-config'   // config-level default
        cache.@repos = 'https://cloud.r-project.org'
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalInstall2rEnv0('dplyr', prefixPath, '--dependencies TRUE')
        then:
        1 * cache.runCommand({ String c ->
            c.contains('--dependencies TRUE') && !c.contains('--from-config') && c.contains('--error')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the install-options override in the env hash' () {
        given:
        def cache = Spy(Install2rCache)
        def BASE = Paths.get('/i2r/envs')

        when:
        def p1 = cache.install2rPrefixPath('dplyr')
        def p2 = cache.install2rPrefixPath('dplyr', '--dependencies TRUE')
        then:
        _ * cache.getCacheDir() >> BASE
        p1 != p2
    }

}
