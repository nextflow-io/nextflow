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

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

class GuixPackageProviderTest extends Specification {

    def 'should report provider name' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        expect:
        provider.getName() == 'guix'
    }

    def 'should only support guix specs' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        expect:
        provider.supportsSpec(new PackageSpec('guix', ['bwa']))
        !provider.supportsSpec(new PackageSpec('conda', ['samtools']))
    }

    def 'should build the activation script' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        when:
        def script = provider.getActivationScript(Paths.get('/work/guix/env-abc'))
        then:
        script.contains('export GUIX_PROFILE=/work/guix/env-abc')
        script.contains('etc/profile')
    }

    def 'should reject a spec for another provider' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('conda', ['x']))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject an empty spec' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('guix', []))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new GuixConfig([:], [:])
        def provider = new GuixPackageProvider(config)
        def cache = Mock(GuixCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('guix', ['bwa', 'samtools']))
        then:
        1 * cache.getCachePathFor('bwa samtools', null) >> Paths.get('/work/guix/env-abc')
    }

    def 'should delegate a single package entry to the cache' () {
        given:
        def config = new GuixConfig([:], [:])
        def provider = new GuixPackageProvider(config)
        def cache = Mock(GuixCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('guix', ['bwa']))
        then:
        1 * cache.getCachePathFor('bwa', null) >> Paths.get('/work/guix/env-abc')
    }

    def 'should delegate an environment-file spec to the cache' () {
        given:
        def config = new GuixConfig([:], [:])
        def provider = new GuixPackageProvider(config)
        def cache = Mock(GuixCache)
        provider.@cache = cache
        // an environment-file spec is a distinct valid config shape from an entries list
        def spec = new PackageSpec('guix').withEnvironment('/path/to/manifest.scm')

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('/path/to/manifest.scm', null) >> Paths.get('/work/guix/env-file')
    }

    def 'should pass per-process install options as an override to the cache' () {
        given:
        def config = new GuixConfig([:], [:])
        def provider = new GuixPackageProvider(config)
        def cache = Mock(GuixCache)
        provider.@cache = cache
        def spec = new PackageSpec('guix', ['bwa'], [installOptions: '--no-grafts'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('bwa', '--no-grafts') >> Paths.get('/work/guix/env-ovr')
    }

    def 'should report config' () {
        given:
        def config = new GuixConfig([:], [:])
        def provider = new GuixPackageProvider(config)
        expect:
        provider.getConfig().is(config)
    }

    def 'should declare manifest file names for auto-detection' () {
        given:
        def provider = new GuixPackageProvider(new GuixConfig([:], [:]))
        expect:
        provider.getManifestFileNames() == ['manifest.scm', 'guix.scm']
    }
}
