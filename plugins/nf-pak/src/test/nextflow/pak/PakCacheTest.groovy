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

package nextflow.pak

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification

class PakCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should create pak env prefix path for a string env' () {
        given:
        def ENV = 'dplyr ggplot2'
        def cache = Spy(PakCache)
        def BASE = Paths.get('/pak/envs')

        when:
        def prefix = cache.pakPrefixPath(ENV)
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/pak/envs/env-')
    }

    def 'should create the correct pak command for package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(PakCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0('dplyr', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('Rscript') && c.contains('pak::pkg_install') && c.contains('"dplyr"')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct pak command for a DESCRIPTION file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def desc = folder.resolve('DESCRIPTION')
        desc.text = 'Package: demo\nImports: dplyr'
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(PakCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0(desc.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('pak::local_install_deps')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should apply a per-process install-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(PakCache)
        cache.@installOptions = 'ask = FALSE'   // config-level default
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0('dplyr', prefixPath, 'upgrade = TRUE')
        then:
        1 * cache.runCommand({ String c ->
            c.contains('upgrade = TRUE') && !c.contains('ask = FALSE')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the install-options override in the env hash' () {
        given:
        def base = Spy(PakCache); base.@installOptions = null
        def ovr  = Spy(PakCache); ovr.@installOptions = null
        def BASE = Paths.get('/pak/envs')

        when:
        def p1 = base.pakPrefixPath('dplyr')
        def p2 = ovr.pakPrefixPath('dplyr', 'upgrade = TRUE')
        then:
        _ * base.getCacheDir() >> BASE
        _ * ovr.getCacheDir() >> BASE
        p1 != p2
    }

    def 'should create the correct pak command for a custom-named manifest file' () {
        // any existing regular file (with a slash in its path) is treated as a
        // manifest and installed via pak::local_install_deps against its parent dir
        given:
        def folder = Files.createTempDirectory('test')
        def manifest = folder.resolve('deps.txt')
        manifest.text = 'Package: demo\nImports: dplyr'
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(PakCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0(manifest.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('pak::local_install_deps') && c.contains(folder.toString())
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should build an R vector with all packages for a multi-package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(PakCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0('dplyr ggplot2 tidyr', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('pak::pkg_install') &&
            c.contains('"dplyr"') && c.contains('"ggplot2"') && c.contains('"tidyr"') &&
            c.contains('c("dplyr", "ggplot2", "tidyr")') &&
            !c.contains('pak::local_install_deps')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'isManifestFile should be true only for an existing regular file with a slash' () {
        given:
        def folder = Files.createTempDirectory('test')
        def file = folder.resolve('deps.txt')
        file.text = 'anything'
        def cache = Spy(PakCache)

        expect:
        // an existing regular file (custom name) is a manifest
        cache.isManifestFile(file.toString())
        // a path with a slash that does not exist is NOT a manifest
        !cache.isManifestFile(folder.resolve('does-not-exist').toString())
        // an existing directory is NOT a regular file, so NOT a manifest
        !cache.isManifestFile(folder.toString())
        // a bare package name (no slash) is NOT a manifest
        !cache.isManifestFile('dplyr')

        cleanup:
        folder?.deleteDir()
    }

    def 'should treat an existing directory spec as a user-provided library directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(PakCache)

        when:
        def prefix = cache.pakPrefixPath(folder.toString())
        then:
        // returned as-is; no cache-dir hashing performed
        0 * cache.getCacheDir()
        prefix == folder

        cleanup:
        folder?.deleteDir()
    }

    def 'should tokenize a multi-line non-file spec into a package list command' () {
        // a multi-line spec that is not an existing file is not a manifest;
        // tokenize() splits on any whitespace including newlines
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(PakCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalPakEnv0('dplyr\nggplot2', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('pak::pkg_install') &&
            c.contains('c("dplyr", "ggplot2")') &&
            !c.contains('pak::local_install_deps')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

}
