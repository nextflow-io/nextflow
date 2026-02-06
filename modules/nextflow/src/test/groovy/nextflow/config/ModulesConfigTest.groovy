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

package nextflow.config

import spock.lang.Specification

/**
 * Tests for ModulesConfig
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModulesConfigTest extends Specification {

    def 'should create empty config'() {
        when:
        def config = new ModulesConfig()

        then:
        config.getAllModules().isEmpty()
    }

    def 'should set and get module version'() {
        given:
        def config = new ModulesConfig()

        when:
        config.setVersion('@nf-core/fastqc', '1.0.0')

        then:
        config.getVersion('@nf-core/fastqc') == '1.0.0'
        config.hasVersion('@nf-core/fastqc')
    }

    def 'should return null for unconfigured module'() {
        given:
        def config = new ModulesConfig()

        when:
        def version = config.getVersion('@nf-core/bwa')

        then:
        version == null
        !config.hasVersion('@nf-core/bwa')
    }

    def 'should override existing version'() {
        given:
        def config = new ModulesConfig()
        config.setVersion('@nf-core/fastqc', '1.0.0')

        when:
        config.setVersion('@nf-core/fastqc', '2.0.0')

        then:
        config.getVersion('@nf-core/fastqc') == '2.0.0'
    }

    def 'should return unmodifiable map from getModules'() {
        given:
        def config = new ModulesConfig()
        config.setVersion('@nf-core/fastqc', '1.0.0')

        when:
        def modules = config.getAllModules()
        modules.put('@nf-core/bwa', '2.0.0')

        then:
        thrown(UnsupportedOperationException)
    }

    def 'should return all configured modules'() {
        given:
        def config = new ModulesConfig()
        config.setVersion('@nf-core/fastqc', '1.0.0')
        config.setVersion('@nf-core/bwa', '2.0.0')
        config.setVersion('@myorg/custom', '0.5.0')

        when:
        def modules = config.getAllModules()

        then:
        modules.size() == 3
        modules['@nf-core/fastqc'] == '1.0.0'
        modules['@nf-core/bwa'] == '2.0.0'
        modules['@myorg/custom'] == '0.5.0'
    }

    def 'should handle empty initialization'() {
        when:
        def config = new ModulesConfig(null)

        then:
        config.getAllModules().isEmpty()
        !config.hasVersion('@nf-core/fastqc')
    }

    def 'should store multiple versions independently'() {
        given:
        def config = new ModulesConfig()

        when:
        config.setVersion('@nf-core/fastqc', '1.0.0')
        config.setVersion('@nf-core/bwa', '2.0.0')
        config.setVersion('@myorg/custom', '0.5.0')

        then:
        config.getVersion('@nf-core/fastqc') == '1.0.0'
        config.getVersion('@nf-core/bwa') == '2.0.0'
        config.getVersion('@myorg/custom') == '0.5.0'
        config.allModules.size() == 3
    }

    def 'should handle module names with special characters'() {
        given:
        def config = new ModulesConfig()

        when:
        config.setVersion('@org-name/module-name', '1.0.0')
        config.setVersion('@org_name/module_name', '2.0.0')
        config.setVersion('simple-module', '3.0.0')

        then:
        config.getVersion('@org-name/module-name') == '1.0.0'
        config.getVersion('@org_name/module_name') == '2.0.0'
        config.getVersion('simple-module') == '3.0.0'
    }

    def 'should handle version strings with various formats'() {
        given:
        def config = new ModulesConfig()

        when:
        config.setVersion('@nf-core/fastqc', '1.0.0')
        config.setVersion('@nf-core/bwa', 'v2.0.0')
        config.setVersion('@nf-core/samtools', '1.0.0-beta')
        config.setVersion('@nf-core/bowtie', '1.0.0-rc.1')

        then:
        config.getVersion('@nf-core/fastqc') == '1.0.0'
        config.getVersion('@nf-core/bwa') == 'v2.0.0'
        config.getVersion('@nf-core/samtools') == '1.0.0-beta'
        config.getVersion('@nf-core/bowtie') == '1.0.0-rc.1'
    }

    def 'should check if multiple modules have versions'() {
        given:
        def config = new ModulesConfig()
        config.setVersion('@nf-core/fastqc', '1.0.0')
        config.setVersion('@nf-core/bwa', '2.0.0')

        expect:
        config.hasVersion('@nf-core/fastqc')
        config.hasVersion('@nf-core/bwa')
        !config.hasVersion('@nf-core/samtools')
    }

    def 'should handle version updates'() {
        given:
        def config = new ModulesConfig()
        config.setVersion('@nf-core/fastqc', '1.0.0')

        when:
        config.setVersion('@nf-core/fastqc', '1.1.0')
        config.setVersion('@nf-core/fastqc', '2.0.0')

        then:
        config.getVersion('@nf-core/fastqc') == '2.0.0'
    }

    def 'should maintain separate versions for different modules'() {
        given:
        def config = new ModulesConfig()

        when:
        config.setVersion('@nf-core/fastqc', '1.0.0')
        config.setVersion('@nf-core/bwa', '2.0.0')

        then:
        config.getVersion('@nf-core/fastqc') == '1.0.0'
        config.getVersion('@nf-core/bwa') == '2.0.0'
        config.getVersion('@nf-core/fastqc') != config.getVersion('@nf-core/bwa')
    }
}
