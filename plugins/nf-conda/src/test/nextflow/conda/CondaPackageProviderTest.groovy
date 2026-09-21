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

package nextflow.conda

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

class CondaPackageProviderTest extends Specification {

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('conda', ['samtools', 'bcftools']))
        then:
        1 * cache.getCachePathFor('samtools bcftools', null) >> Paths.get('/work/conda/env-abc')
    }

    def 'should pass per-process create options as an override to the cache' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache
        def spec = new PackageSpec('conda', ['samtools'], [createOptions: '--override-channels'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('samtools', '--override-channels') >> Paths.get('/work/conda/env-ovr')
    }

    def 'should delegate a single package entry to the cache' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('conda', ['samtools=1.17']))
        then:
        1 * cache.getCachePathFor('samtools=1.17', null) >> Paths.get('/work/conda/env-one')
    }

    def 'should delegate an environment file to the cache' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache
        def spec = new PackageSpec('conda').withEnvironment('/some/environment.yml')

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('/some/environment.yml', null) >> Paths.get('/work/conda/env-file')
    }

    def 'should throw for an empty spec with neither entries nor environment' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('conda'))
        then:
        thrown(IllegalArgumentException)
        0 * cache.getCachePathFor(_, _)
    }

    def 'should throw when the spec targets a different provider' () {
        given:
        def config = new CondaConfig([:], [:])
        def provider = new CondaPackageProvider(config)
        def cache = Mock(CondaCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('pixi', ['samtools']))
        then:
        thrown(IllegalArgumentException)
        0 * cache.getCachePathFor(_, _)
    }
}
