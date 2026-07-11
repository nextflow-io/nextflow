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

package nextflow.pixi

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

class PixiPackageProviderTest extends Specification {

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('pixi', ['cowpy', 'numpy']))
        then:
        1 * cache.getCachePathFor('cowpy numpy', null) >> Paths.get('/work/pixi/env-abc')
    }

    def 'should pass per-process create options as an override to the cache' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache
        def spec = new PackageSpec('pixi', ['cowpy'], [createOptions: '--frozen'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('cowpy', '--frozen') >> Paths.get('/work/pixi/env-ovr')
    }

    def 'should delegate a single-package entry to the cache' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('pixi', ['cowpy']))
        then:
        1 * cache.getCachePathFor('cowpy', null) >> Paths.get('/work/pixi/env-single')
    }

    def 'should delegate an environment-file spec to the cache' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache
        def spec = new PackageSpec('pixi').withEnvironment('/path/to/myenv.toml')

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('/path/to/myenv.toml', null) >> Paths.get('/work/pixi/env-file')
    }

    def 'should throw for an empty spec with no entries and no environment' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('pixi'))
        then:
        thrown(IllegalArgumentException)
        and:
        0 * cache.getCachePathFor(_, _)
    }

    def 'should throw when the spec provider is a different provider' () {
        given:
        def config = new PixiConfig([:], [:])
        def provider = new PixiPackageProvider(config)
        def cache = Mock(PixiCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('conda', ['cowpy']))
        then:
        thrown(IllegalArgumentException)
        and:
        0 * cache.getCachePathFor(_, _)
    }
}
