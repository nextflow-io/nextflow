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

import java.nio.file.Files

import nextflow.ISession
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

    def 'should parse the directive named-args form (leading options map)'() {
        when:
        // `package "numpy pandas", provider: "uv"` is delivered as [ [provider:'uv'], "numpy pandas" ]
        def spec = PackageManager.parseSpec([[provider: 'uv'], 'numpy pandas'], 'conda')

        then:
        spec.provider == 'uv'
        spec.entries == ['numpy pandas']
    }

    def 'should fall back to the default provider when none is given in the named-args form'() {
        when:
        def spec = PackageManager.parseSpec([[channels: ['bioconda']], 'samtools'], 'conda')

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

    def 'should auto-detect a manifest file in the module directory'() {
        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('environment.yml').text = 'name: test'
        def manifests = [conda: ['environment.yml', 'environment.yaml'], uv: ['requirements.txt']]

        when:
        def spec = PackageManager.findManifestSpec(dir, manifests)

        then:
        spec.provider == 'conda'
        spec.environment == dir.resolve('environment.yml').toAbsolutePath().toString()
        spec.hasEnvironmentFile()

        cleanup:
        dir?.deleteDir()
    }

    def 'should auto-detect the uv manifest when only it is present'() {
        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('requirements.txt').text = 'numpy'
        def manifests = [conda: ['environment.yml'], uv: ['requirements.txt', 'pyproject.toml']]

        when:
        def spec = PackageManager.findManifestSpec(dir, manifests)

        then:
        spec.provider == 'uv'
        spec.environment.endsWith('requirements.txt')

        cleanup:
        dir?.deleteDir()
    }

    def 'should prefer conda over uv when both manifests are present'() {
        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('environment.yml').text = 'name: test'
        dir.resolve('requirements.txt').text = 'numpy'
        def manifests = [conda: ['environment.yml'], uv: ['requirements.txt']]

        when:
        def spec = PackageManager.findManifestSpec(dir, manifests)

        then:
        spec.provider == 'conda'

        cleanup:
        dir?.deleteDir()
    }

    def 'should return null when no manifest is found'() {
        given:
        def dir = Files.createTempDirectory('test')
        def manifests = [conda: ['environment.yml'], uv: ['requirements.txt']]

        expect:
        PackageManager.findManifestSpec(dir, manifests) == null

        cleanup:
        dir?.deleteDir()
    }

    def 'should return null when module dir is null'() {
        expect:
        PackageManager.findManifestSpec(null, [conda: ['environment.yml']]) == null
    }

    def 'should check if the feature is enabled'() {
        given:
        // use a real config map so the `navigate` extension method resolves
        // (mocking Map fails because `navigate` is a Groovy extension method)
        def session = Mock(ISession) { getConfig() >> config }

        expect:
        PackageManager.isEnabled(session) == result

        where:
        config                                        | result
        [nextflow: [preview: [package: true]]]        | true
        [nextflow: [preview: [package: false]]]       | false
        [nextflow: [preview: [:]]]                    | false
        [:]                                           | false
    }

    def 'should check if the feature is enabled by the script feature flag'() {
        given:
        def session = Mock(ISession) { getConfig() >> [:] }

        when: 'the `nextflow.preview.package` flag is set in the pipeline script'
        nextflow.NextflowMeta.instance.preview.setPackage(true)
        then:
        PackageManager.isEnabled(session)

        cleanup:
        nextflow.NextflowMeta.instance.preview.setPackage(false)
    }
}