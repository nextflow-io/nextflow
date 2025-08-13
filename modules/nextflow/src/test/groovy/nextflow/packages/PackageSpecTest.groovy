/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.packages

import spock.lang.Specification

/**
 * Unit tests for PackageSpec
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class PackageSpecTest extends Specification {

    def 'should create package spec with builder pattern'() {
        when:
        def spec = new PackageSpec()
            .withProvider('conda')
            .withEntries(['samtools=1.17', 'bcftools=1.18'])
            .withChannels(['conda-forge', 'bioconda'])

        then:
        spec.provider == 'conda'
        spec.entries == ['samtools=1.17', 'bcftools=1.18']
        spec.channels == ['conda-forge', 'bioconda']
        spec.isValid()
        spec.hasEntries()
        !spec.hasEnvironmentFile()
    }

    def 'should create package spec with environment file'() {
        when:
        def spec = new PackageSpec()
            .withProvider('conda')
            .withEnvironment('name: myenv\ndependencies:\n  - samtools=1.17')

        then:
        spec.provider == 'conda'
        spec.environment == 'name: myenv\ndependencies:\n  - samtools=1.17'
        spec.isValid()
        !spec.hasEntries()
        spec.hasEnvironmentFile()
    }

    def 'should validate spec correctly'() {
        expect:
        new PackageSpec().withProvider('conda').withEntries(['samtools']).isValid()
        new PackageSpec().withProvider('conda').withEnvironment('deps').isValid()
        !new PackageSpec().withProvider('conda').isValid() // no entries or environment
        !new PackageSpec().withEntries(['samtools']).isValid() // no provider
    }

    def 'should handle constructor with parameters'() {
        when:
        def spec = new PackageSpec('pixi', ['samtools=1.17'], [channels: ['conda-forge']])

        then:
        spec.provider == 'pixi'
        spec.entries == ['samtools=1.17']
        spec.options == [channels: ['conda-forge']]
    }

    def 'should handle empty or null values'() {
        when:
        def spec = new PackageSpec('conda', null, null)

        then:
        spec.provider == 'conda'
        spec.entries == []
        spec.options == [:]
    }
}