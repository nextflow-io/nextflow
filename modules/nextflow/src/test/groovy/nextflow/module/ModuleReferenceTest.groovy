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

package nextflow.module

import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 * Test suite for ModuleReference
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleReferenceTest extends Specification {

    def 'should parse valid module reference with @'() {
        when:
        def ref = ModuleReference.parse('@nf-core/fastqc')

        then:
        ref.scope == 'nf-core'
        ref.name == 'fastqc'
        ref.fullName == '@nf-core/fastqc'
    }

    def 'should parse valid module reference without @'() {
        when:
        def ref = ModuleReference.parse('nf-core/fastqc')

        then:
        ref.scope == 'nf-core'
        ref.name == 'fastqc'
        ref.fullName == '@nf-core/fastqc'
    }

    def 'should parse module reference with multiple slashes'() {
        when:
        def ref = ModuleReference.parse('@myorg/samtools/view')

        then:
        ref.scope == 'myorg'
        ref.name == 'samtools/view'
        ref.fullName == '@myorg/samtools/view'
    }

    def 'should reject invalid module reference without scope'() {
        when:
        ModuleReference.parse('fastqc')

        then:
        thrown(AbortOperationException)
    }

    def 'should reject empty module reference'() {
        when:
        ModuleReference.parse('')

        then:
        thrown(AbortOperationException)
    }

    def 'should reject null module reference'() {
        when:
        ModuleReference.parse(null)

        then:
        thrown(AbortOperationException)
    }

    def 'should reject module reference with only @'() {
        when:
        ModuleReference.parse('@')

        then:
        thrown(AbortOperationException)
    }

    def 'should reject module reference with only scope'() {
        when:
        ModuleReference.parse('@nf-core/')

        then:
        thrown(AbortOperationException)
    }

    def 'should handle module reference with trailing slash'() {
        when:
        ModuleReference.parse('@nf-core/fastqc/')

        then:
        thrown(AbortOperationException)
    }

    def 'should create module reference from components'() {
        when:
        def ref = new ModuleReference('nf-core', 'fastqc')

        then:
        ref.scope == 'nf-core'
        ref.name == 'fastqc'
        ref.fullName == '@nf-core/fastqc'
    }

    def 'should handle scope names with hyphens'() {
        when:
        def ref = ModuleReference.parse('@my-org/my-module')

        then:
        ref.scope == 'my-org'
        ref.name == 'my-module'
    }

    def 'should handle scope names with underscores'() {
        when:
        def ref = ModuleReference.parse('@my_org/my_module')

        then:
        ref.scope == 'my_org'
        ref.name == 'my_module'
    }

    def 'should handle module names with numbers'() {
        when:
        def ref = ModuleReference.parse('@nf-core/bwa-mem2')

        then:
        ref.scope == 'nf-core'
        ref.name == 'bwa-mem2'
    }

    def 'should implement equals correctly'() {
        given:
        def ref1 = ModuleReference.parse('@nf-core/fastqc')
        def ref2 = ModuleReference.parse('@nf-core/fastqc')
        def ref3 = ModuleReference.parse('@nf-core/multiqc')

        expect:
        ref1 == ref2
        ref1 != ref3
    }

    def 'should implement hashCode correctly'() {
        given:
        def ref1 = ModuleReference.parse('@nf-core/fastqc')
        def ref2 = ModuleReference.parse('@nf-core/fastqc')

        expect:
        ref1.hashCode() == ref2.hashCode()
    }

    def 'should implement toString correctly'() {
        given:
        def ref = ModuleReference.parse('@nf-core/fastqc')

        expect:
        ref.toString() == '@nf-core/fastqc'
    }

    def 'should be usable as map key'() {
        given:
        def ref1 = ModuleReference.parse('@nf-core/fastqc')
        def ref2 = ModuleReference.parse('@nf-core/fastqc')
        def ref3 = ModuleReference.parse('@nf-core/multiqc')

        def map = [:]
        map[ref1] = 'value1'
        map[ref3] = 'value3'

        expect:
        map[ref2] == 'value1'  // ref2 equals ref1, should get same value
        map[ref3] == 'value3'
        map.size() == 2
    }

    def 'should handle org-style scopes'() {
        when:
        def ref = ModuleReference.parse('@mycompany.io/custom-module')

        then:
        ref.scope == 'mycompany.io'
        ref.name == 'custom-module'
    }

    def 'should reject module reference with spaces'() {
        when:
        ModuleReference.parse('@nf-core/fast qc')

        then:
        thrown(AbortOperationException)
    }

    def 'should reject module reference with special characters'() {
        when:
        ModuleReference.parse('@nf-core/fastqc!')

        then:
        thrown(AbortOperationException)
    }

    def 'should handle deeply nested module names'() {
        when:
        def ref = ModuleReference.parse('@nf-core/samtools/sort/parallel')

        then:
        ref.scope == 'nf-core'
        ref.name == 'samtools/sort/parallel'
        ref.fullName == '@nf-core/samtools/sort/parallel'
    }

    def 'should parse from string with leading/trailing whitespace'() {
        when:
        def ref = ModuleReference.parse('  @nf-core/fastqc  ')

        then:
        ref.scope == 'nf-core'
        ref.name == 'fastqc'
    }
}
