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

package nextflow.nix

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification

class NixCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should create nix env prefix path for a string env' () {
        given:
        def ENV = 'bwa samtools'
        def cache = Spy(NixCache)
        def BASE = Paths.get('/nix/envs')

        when:
        def prefix = cache.nixPrefixPath(ENV)
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/nix/envs/env-')
    }

    def 'should include flake ref in hash' () {
        given:
        def cache1 = Spy(NixCache)
        cache1.@flakeRef = 'nixpkgs'
        def cache2 = Spy(NixCache)
        cache2.@flakeRef = 'other'
        def BASE = Paths.get('/nix/envs')

        when:
        def prefix1 = cache1.nixPrefixPath('bwa')
        def prefix2 = cache2.nixPrefixPath('bwa')
        then:
        _ * cache1.getCacheDir() >> BASE
        _ * cache2.getCacheDir() >> BASE
        prefix1 != prefix2
    }

    def 'should create the correct nix command for package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalNixEnv0('bwa', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('nix') && c.contains('profile install') && c.contains('nixpkgs#bwa')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should build the nix command for a multiple-package list' () {
        // NOTE: nix does NOT support manifest files. A package list is space-separated
        // and each bare token is mapped to `<flakeRef>#<token>`; there is no custom-named
        // manifest-file success path to test (a slash-path is treated as an installable,
        // see the 'path-with-slash' test below).
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-multi')
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalNixEnv0('hello samtools', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('profile install') && c.contains('nixpkgs#hello') && c.contains('nixpkgs#samtools')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should not prepend the flake ref to a token that already has one' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-flakeref')
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        // a fully-qualified installable already contains '#' and must pass through unchanged
        cache.createLocalNixEnv0('github:owner/repo#hello', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('github:owner/repo#hello') && !c.contains('nixpkgs#github:')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should hash a slash-path that is not an existing directory as an installable, not a manifest file' () {
        // nix has no manifest-file concept: a path-like spec that is not an existing
        // directory falls through to the normal installable-tokenisation path and is
        // hashed into a cache dir (it is NOT read as a file).
        given:
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        def BASE = Paths.get('/nix/envs')

        when:
        def prefix = cache.nixPrefixPath('/does/not/exist/env.nix')
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/nix/envs/env-')
    }

    def 'should pass a slash-path straight through as an installable in the command' () {
        // the path is tokenised and prefixed with the flake ref (NOT read as a manifest)
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-slash')
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalNixEnv0('/does/not/exist/env.nix', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('nixpkgs#/does/not/exist/env.nix')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should return an existing profile directory as-is without hashing' () {
        given:
        def folder = Files.createTempDirectory('test')
        def profileDir = folder.resolve('my-profile')
        Files.createDirectory(profileDir)
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'

        when:
        def prefix = cache.nixPrefixPath(profileDir.toString())
        then:
        // an existing directory path is used directly as the profile, never cache-hashed
        0 * cache.getCacheDir()
        prefix == profileDir

        cleanup:
        folder?.deleteDir()
    }

    def 'should apply a per-process install-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(NixCache)
        cache.@flakeRef = 'nixpkgs'
        cache.@installOptions = '--from-config'   // config-level default
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalNixEnv0('bwa', prefixPath, '--offline')
        then:
        1 * cache.runCommand({ String c ->
            c.contains('--offline') && !c.contains('--from-config')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the install-options override in the env hash' () {
        given:
        def base = Spy(NixCache); base.@flakeRef = 'nixpkgs'
        def ovr  = Spy(NixCache); ovr.@flakeRef = 'nixpkgs'
        def BASE = Paths.get('/nix/envs')

        when:
        def p1 = base.nixPrefixPath('bwa')
        def p2 = ovr.nixPrefixPath('bwa', '--offline')
        then:
        _ * base.getCacheDir() >> BASE
        _ * ovr.getCacheDir() >> BASE
        p1 != p2
    }

}
