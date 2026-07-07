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

import java.nio.file.Paths

import nextflow.packages.PackageSpec
import spock.lang.Specification

/**
 * NOTE: install2.r has no concept of a manifest/environment file -- packages are
 * always passed as bare names on the command line. There is therefore no
 * custom-named-manifest success test here; a path is only ever treated as a
 * package name (or, in the cache, as a verbatim library dir when it exists).
 */
class Install2rPackageProviderTest extends Specification {

    def 'should report provider name' () {
        given:
        def provider = new Install2rPackageProvider(new Install2rConfig([:], [:]))
        expect:
        provider.getName() == 'install2r'
    }

    def 'should only support install2r specs' () {
        given:
        def provider = new Install2rPackageProvider(new Install2rConfig([:], [:]))
        expect:
        provider.supportsSpec(new PackageSpec('install2r', ['dplyr']))
        !provider.supportsSpec(new PackageSpec('conda', ['samtools']))
    }

    def 'should build the activation script' () {
        given:
        def provider = new Install2rPackageProvider(new Install2rConfig([:], [:]))
        when:
        def script = provider.getActivationScript(Paths.get('/work/install2r/env-abc'))
        then:
        script.contains('export R_LIBS_USER=/work/install2r/env-abc')
    }

    def 'should reject a spec for another provider' () {
        given:
        def provider = new Install2rPackageProvider(new Install2rConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('conda', ['samtools']))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject an empty spec' () {
        given:
        def provider = new Install2rPackageProvider(new Install2rConfig([:], [:]))
        when:
        provider.createEnvironment(new PackageSpec('install2r', []))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should delegate package entries to the cache as a space-separated env' () {
        given:
        def config = new Install2rConfig([:], [:])
        def provider = new Install2rPackageProvider(config)
        def cache = Mock(Install2rCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('install2r', ['dplyr', 'ggplot2']))
        then:
        1 * cache.getCachePathFor('dplyr ggplot2', null) >> Paths.get('/work/install2r/env-abc')
    }

    def 'should delegate a multi-package list joined by spaces to the cache' () {
        given:
        def config = new Install2rConfig([:], [:])
        def provider = new Install2rPackageProvider(config)
        def cache = Mock(Install2rCache)
        provider.@cache = cache

        when:
        provider.createEnvironment(new PackageSpec('install2r', ['dplyr', 'ggplot2', 'tidyr']))
        then:
        1 * cache.getCachePathFor('dplyr ggplot2 tidyr', null) >> Paths.get('/work/install2r/env-multi')
    }

    def 'should pass per-process install options as an override to the cache' () {
        given:
        def config = new Install2rConfig([:], [:])
        def provider = new Install2rPackageProvider(config)
        def cache = Mock(Install2rCache)
        provider.@cache = cache
        def spec = new PackageSpec('install2r', ['dplyr'], [installOptions: '--repos http://mirror'])

        when:
        provider.createEnvironment(spec)
        then:
        1 * cache.getCachePathFor('dplyr', '--repos http://mirror') >> Paths.get('/work/install2r/env-ovr')
    }

    def 'should report config' () {
        given:
        def config = new Install2rConfig([:], [:])
        def provider = new Install2rPackageProvider(config)
        expect:
        provider.getConfig().is(config)
    }
}
