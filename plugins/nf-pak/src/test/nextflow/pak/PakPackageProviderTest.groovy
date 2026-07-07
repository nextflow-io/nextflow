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

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

class PakPackageProviderTest extends Specification {

    def 'should report provider name' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        expect:
        provider.getName() == 'pak'
    }

    def 'should only support pak specs' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        expect:
        provider.supportsSpec(new PackageSpec('pak', ['dplyr']))
        !provider.supportsSpec(new PackageSpec('conda', ['samtools']))
    }

    def 'should build the activation script' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        when:
        def script = provider.getActivationScript(Paths.get('/work/pak/env-abc'))
        then:
        script.contains('export R_LIBS_USER=/work/pak/env-abc')
    }

    def 'should reject a spec for another provider' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('conda', ['samtools']))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject an empty spec' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('pak', []))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new PakConfig([:], [:])
        def provider = new PakPackageProvider(config)
        def cache = Mock(PakCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('pak', ['dplyr', 'ggplot2']))
        then:
        1 * cache.getCachePathFor('dplyr ggplot2', null) >> Paths.get('/work/pak/env-abc')
    }

    def 'should delegate an environment library directory to the cache' () {
        given:
        def config = new PakConfig([:], [:])
        def provider = new PakPackageProvider(config)
        def cache = Mock(PakCache)
        provider.@cache = cache
        def spec = new PackageSpec('pak', []).withEnvironment('/opt/R/lib')

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('/opt/R/lib', null) >> Paths.get('/opt/R/lib')
    }

    def 'should pass per-process install options as an override to the cache' () {
        given:
        def config = new PakConfig([:], [:])
        def provider = new PakPackageProvider(config)
        def cache = Mock(PakCache)
        provider.@cache = cache
        def spec = new PackageSpec('pak', ['dplyr'], [installOptions: 'upgrade = TRUE'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('dplyr', 'upgrade = TRUE') >> Paths.get('/work/pak/env-ovr')
    }

    def 'should report config' () {
        given:
        def config = new PakConfig([:], [:])
        def provider = new PakPackageProvider(config)
        expect:
        provider.getConfig().is(config)
    }

    def 'should declare manifest file names for auto-detection' () {
        given:
        def provider = new PakPackageProvider(new PakConfig([:], [:]))
        expect:
        provider.getManifestFileNames() == ['DESCRIPTION']
    }
}
