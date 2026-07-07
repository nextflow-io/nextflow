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

package nextflow.guix

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import nextflow.util.Duration
import spock.lang.Specification

class GuixCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should create guix env prefix path for a string env' () {
        given:
        def ENV = 'bwa samtools'
        def cache = Spy(GuixCache)
        def BASE = Paths.get('/guix/envs')

        when:
        def prefix = cache.guixPrefixPath(ENV)
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/guix/envs/env-')
    }

    def 'should create the correct guix command for package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalGuixEnv0('bwa samtools', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--install') && c.contains('bwa samtools')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should detect a manifest file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def manifest = folder.resolve('manifest.scm')
        manifest.text = '(specifications->manifest (list "bwa"))'
        def cache = Spy(GuixCache)

        expect:
        cache.isManifestFile(manifest.toString())
        !cache.isManifestFile('bwa samtools')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct guix command for a manifest file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def manifest = folder.resolve('manifest.scm')
        manifest.text = '(specifications->manifest (list "bwa"))'
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalGuixEnv0(manifest.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--manifest=') && c.contains('manifest.scm')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct guix command for a custom-named manifest file' () {
        given:
        def folder = Files.createTempDirectory('test')
        // a NON-conventional name with the supported extension proves the name is not special
        def manifest = folder.resolve('deps.scm')
        manifest.text = '(specifications->manifest (list "bwa"))'
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalGuixEnv0(manifest.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--manifest=') && c.contains('deps.scm') && !c.contains('--install')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should detect a custom-named regular file as a manifest' () {
        given:
        def folder = Files.createTempDirectory('test')
        def manifest = folder.resolve('deps.scm')
        manifest.text = '(specifications->manifest (list "bwa"))'
        def cache = Spy(GuixCache)

        expect:
        // isManifestFile is name-agnostic: any existing regular file at a path-with-slash
        cache.isManifestFile(manifest.toString())

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct guix command for a multi-package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-multi')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalGuixEnv0('bwa samtools bcftools', prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--install') && c.contains('bwa samtools bcftools') && !c.contains('--manifest=')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should treat a path-with-slash that is not an existing directory as a package spec' () {
        given:
        def cache = Spy(GuixCache)
        def BASE = Paths.get('/guix/envs')
        // a slash-containing spec that does not resolve to an existing directory or file
        def spec = '/no/such/dir/bwa'

        when:
        def prefix = cache.guixPrefixPath(spec)
        then:
        1 * cache.getCacheDir() >> BASE
        // not a manifest file, so it is hashed as an env string and lands under the cache dir
        !cache.isManifestFile(spec)
        prefix.toString().startsWith('/guix/envs/env-')
    }

    def 'should create the correct guix command for a path-with-slash that is not a directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')
        def spec = '/no/such/dir/bwa'

        when:
        cache.createLocalGuixEnv0(spec, prefixPath)
        then:
        // not a regular file -> falls through to --install rather than --manifest
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--install') && c.contains(spec) && !c.contains('--manifest=')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct guix command for a multi-line spec' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-x')
        def cache = Spy(GuixCache)
        cache.@installOptions = null
        cache.@createTimeout = Duration.of('20min')
        // a multi-line spec has no slash and is not a file -> treated as an --install string
        def spec = 'bwa\nsamtools'

        when:
        cache.createLocalGuixEnv0(spec, prefixPath)
        then:
        1 * cache.runCommand({ String c ->
            c.contains('guix package') && c.contains('--install') && !c.contains('--manifest=')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should apply a per-process install-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(GuixCache)
        cache.@installOptions = '--from-config'   // config-level default
        cache.@createTimeout = Duration.of('20min')

        when:
        cache.createLocalGuixEnv0('bwa', prefixPath, '--no-grafts')
        then:
        1 * cache.runCommand({ String c ->
            c.contains('--no-grafts') && !c.contains('--from-config')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the install-options override in the env hash' () {
        given:
        def base = Spy(GuixCache); base.@installOptions = null
        def ovr  = Spy(GuixCache); ovr.@installOptions = null
        def BASE = Paths.get('/guix/envs')

        when:
        def p1 = base.guixPrefixPath('bwa')
        def p2 = ovr.guixPrefixPath('bwa', '--no-grafts')
        then:
        _ * base.getCacheDir() >> BASE
        _ * ovr.getCacheDir() >> BASE
        p1 != p2
    }

    def 'should hash manifest content rather than the path' () {
        given:
        def folder = Files.createTempDirectory('test')
        def manifest = folder.resolve('manifest.scm')
        manifest.text = '(specifications->manifest (list "bwa"))'
        def cache = Spy(GuixCache)
        def BASE = Paths.get('/guix/envs')

        when:
        def prefix = cache.guixPrefixPath(manifest.toString())
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/guix/envs/env-')

        cleanup:
        folder?.deleteDir()
    }

}
