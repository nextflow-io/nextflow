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

import nextflow.Session
import spock.lang.Specification

/**
 * Unit tests for PackageManager
 * 
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class PackageManagerTest extends Specification {

    def 'should parse string package definition'() {
        when:
        def spec = PackageManager.parseSpec('samtools=1.17', 'conda')

        then:
        spec.provider == 'conda'
        spec.entries == ['samtools=1.17']
    }

    def 'should parse list package definition'() {
        when:
        def spec = PackageManager.parseSpec(['samtools=1.17', 'bcftools=1.18'], 'conda')

        then:
        spec.provider == 'conda'
        spec.entries == ['samtools=1.17', 'bcftools=1.18']
    }

    def 'should parse map package definition'() {
        when:
        def spec = PackageManager.parseSpec([
            provider: 'pixi',
            packages: ['samtools=1.17', 'bcftools=1.18'],
            channels: ['conda-forge', 'bioconda']
        ], 'conda')

        then:
        spec.provider == 'pixi'
        spec.entries == ['samtools=1.17', 'bcftools=1.18']
        spec.channels == ['conda-forge', 'bioconda']
    }

    def 'should parse map with environment file'() {
        when:
        def spec = PackageManager.parseSpec([
            provider: 'conda',
            environment: 'name: test\ndependencies:\n  - samtools'
        ], null)

        then:
        spec.provider == 'conda'
        spec.environment == 'name: test\ndependencies:\n  - samtools'
        spec.hasEnvironmentFile()
    }

    def 'should use default provider when not specified'() {
        when:
        def spec = PackageManager.parseSpec([
            packages: ['samtools=1.17']
        ], 'conda')

        then:
        spec.provider == 'conda'
        spec.entries == ['samtools=1.17']
    }

    def 'should handle single package in map'() {
        when:
        def spec = PackageManager.parseSpec([
            provider: 'pixi',
            packages: 'samtools=1.17'
        ], null)

        then:
        spec.provider == 'pixi'
        spec.entries == ['samtools=1.17']
    }

    def 'should handle single channel in map'() {
        when:
        def spec = PackageManager.parseSpec([
            provider: 'conda',
            packages: ['samtools'],
            channels: 'bioconda'
        ], null)

        then:
        spec.provider == 'conda'
        spec.entries == ['samtools']
        spec.channels == ['bioconda']
    }

    def 'should throw error for invalid package definition'() {
        when:
        PackageManager.parseSpec(123, 'conda')

        then:
        thrown(IllegalArgumentException)
    }

    def 'should check if feature is enabled'() {
        given:
        def session = Mock(Session) {
            getConfig() >> Mock() {
                navigate('nextflow.preview.package', false) >> enabled
            }
        }

        expect:
        PackageManager.isEnabled(session) == result

        where:
        enabled | result
        true    | true
        false   | false
        null    | false
    }
}