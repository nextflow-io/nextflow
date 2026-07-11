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

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

class NixPackageProviderTest extends Specification {

    def 'should report provider name' () {
        given:
        def provider = new NixPackageProvider(new NixConfig([:], [:]))
        expect:
        provider.getName() == 'nix'
    }

    def 'should only support nix specs' () {
        given:
        def provider = new NixPackageProvider(new NixConfig([:], [:]))
        expect:
        provider.supportsSpec(new PackageSpec('nix', ['bwa']))
        !provider.supportsSpec(new PackageSpec('conda', ['x']))
    }

    def 'should build the activation script' () {
        given:
        def provider = new NixPackageProvider(new NixConfig([:], [:]))
        when:
        def script = provider.getActivationScript(Paths.get('/work/nix/env-abc'))
        then:
        script.contains('export PATH=/work/nix/env-abc/bin')
    }

    def 'should reject a spec for another provider' () {
        given:
        def provider = new NixPackageProvider(new NixConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('conda', ['x']))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject an empty spec' () {
        given:
        def provider = new NixPackageProvider(new NixConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('nix', []))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new NixConfig([:], [:])
        def provider = new NixPackageProvider(config)
        def cache = Mock(NixCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('nix', ['bwa', 'samtools']))
        then:
        1 * cache.getCachePathFor('bwa samtools', null) >> Paths.get('/work/nix/env-abc')
    }

    def 'should delegate a single-package entry to the cache' () {
        // second valid-config shape: a single-element package list
        given:
        def config = new NixConfig([:], [:])
        def provider = new NixPackageProvider(config)
        def cache = Mock(NixCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('nix', ['hello']))
        then:
        1 * cache.getCachePathFor('hello', null) >> Paths.get('/work/nix/env-hello')
    }

    def 'should delegate an environment-file (profile directory) spec to the cache' () {
        // valid-config shape: an existing Nix profile directory supplied via `environment`.
        // NOTE: nix does NOT support manifest files -- `environment` here is a profile dir
        // path, which is forwarded verbatim to the cache (the cache decides whether it is an
        // existing directory or an installable).
        given:
        def config = new NixConfig([:], [:])
        def provider = new NixPackageProvider(config)
        def cache = Mock(NixCache)
        provider.@cache = cache
        def spec = new PackageSpec('nix').withEnvironment('/opt/nix/my-profile')

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('/opt/nix/my-profile', null) >> Paths.get('/opt/nix/my-profile')
    }

    def 'should pass per-process install options as an override to the cache' () {
        given:
        def config = new NixConfig([:], [:])
        def provider = new NixPackageProvider(config)
        def cache = Mock(NixCache)
        provider.@cache = cache
        def spec = new PackageSpec('nix', ['bwa'], [installOptions: '--offline'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('bwa', '--offline') >> Paths.get('/work/nix/env-ovr')
    }

    def 'should report config' () {
        given:
        def config = new NixConfig([:], [:])
        def provider = new NixPackageProvider(config)
        expect:
        provider.getConfig().is(config)
    }
}
