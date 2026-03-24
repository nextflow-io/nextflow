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
import spock.lang.TempDir

import java.nio.file.Path

/**
 * Tests for ModuleSpec
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleSpecTest extends Specification {

    @TempDir
    Path tempDir

    def 'should load valid spec' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/fastqc
version: 1.0.0
description: FastQC quality control
authors:
  - John Doe
license: MIT
keywords:
  - quality-control
  - fastq
requires:
  nextflow: ">=24.04.0"
'''

        when:
        def spec = ModuleSpec.load(metaYaml)

        then:
        spec.name == 'nf-core/fastqc'
        spec.version == '1.0.0'
        spec.description == 'FastQC quality control'
        spec.authors == ['John Doe']
        spec.license == 'MIT'
        spec.keywords == ['quality-control', 'fastq']
        spec.requires == ['nextflow': '>=24.04.0']
    }

    def 'should fail to load non-existent spec' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')

        when:
        ModuleSpec.load(metaYaml)

        then:
        thrown(AbortOperationException)
    }

    def 'should validate complete spec' () {
        given:
        def spec = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: '1.0.0',
            description: 'FastQC quality control',
            license: 'MIT'
        )

        when:
        def errors = spec.validate()

        then:
        errors.isEmpty()
        spec.isValid()
    }

    def 'should detect missing required fields' () {
        given:
        def spec = new ModuleSpec(
            name: 'nf-core/fastqc'
            // missing version, description, license
        )

        when:
        def errors = spec.validate()

        then:
        errors.size() == 3
        errors.any { it.contains('version') }
        errors.any { it.contains('description') }
        errors.any { it.contains('license') }
        !spec.isValid()
    }

    def 'should validate version format' () {
        given:
        def spec = new ModuleSpec(
            name: 'nf-core/fastqc',
            version: version,
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = spec.validate()

        then:
        errors.isEmpty() == valid

        where:
        version         | valid
        '1.0.0'         | true
        '1.0.0-alpha'   | true
        '1.0.0-beta.1'  | true
        '1.0'           | false
        'v1.0.0'        | false
        '1.0.0.0'       | false
    }

    def 'should return null from loadInputTypes when file does not exist' () {
        expect:
        ModuleSpec.loadInputTypes(tempDir.resolve('nonexistent.yml')) == null
    }

    def 'should return null from loadInputTypes when no input section' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/fastqc
description: FastQC quality control
'''
        expect:
        ModuleSpec.loadInputTypes(metaYaml) == null
    }

    def 'should load input types in new paramSpec format' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/fastqc
description: FastQC quality control
input:
  - name: meta
    type: map
    description: Sample metadata
  - name: reads
    type: file
    description: Input reads
'''
        when:
        def types = ModuleSpec.loadInputTypes(metaYaml)

        then:
        types == [meta: 'map', reads: 'file']
    }

    def 'should load input types in old nf-core format' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/fastqc
description: FastQC quality control
input:
  - - meta:
        type: map
        description: Sample metadata
    - reads:
        type: file
        description: Input reads
  - index:
      type: directory
      description: Index directory
'''
        when:
        def types = ModuleSpec.loadInputTypes(metaYaml)

        then:
        types == [meta: 'map', reads: 'file', index: 'directory']
    }

    def 'should flatten tuple inputs in paramSpec format' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/tool
description: Tool
input:
  - - name: meta
      type: map
      description: Sample metadata
    - name: reads
      type: file
      description: Input reads
  - name: index
    type: directory
    description: Index directory
'''
        when:
        def types = ModuleSpec.loadInputTypes(metaYaml)

        then:
        types == [meta: 'map', reads: 'file', index: 'directory']
    }

    def 'should flatten tuple inputs in nf-core tuple-as-map-value format' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = '''
name: nf-core/tool
description: Tool
input:
  - tuple:
      - meta:
          type: map
          description: Sample metadata
      - reads:
          type: file
          description: Input reads
  - index:
      type: directory
      description: Index directory
'''
        when:
        def types = ModuleSpec.loadInputTypes(metaYaml)

        then:
        types == [meta: 'map', reads: 'file', index: 'directory']
    }

    def 'should return null from loadInputTypes for malformed YAML' () {
        given:
        def metaYaml = tempDir.resolve('meta.yml')
        metaYaml.text = ': invalid: yaml: {'

        expect:
        ModuleSpec.loadInputTypes(metaYaml) == null
    }

    def 'should validate module name format' () {
        given:
        def spec = new ModuleSpec(
            name: name,
            version: '1.0.0',
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = spec.validate()

        then:
        errors.isEmpty() == valid

        where:
        name                        | valid
        'nf-core/fastqc'            | true
        'myorg/my-module'           | true
        'org_1/tool_2'              | true
        'nf-core/gfatools/gfa2fa'   | true   // nested module path
        'myorg/tools/sub/module'    | true   // deeply nested
        'org.name/tool/sub'         | true   // dot in scope
        'fastqc'                    | false
        '@nf-core/fastqc'           | false
        'nf-core/fast qc'           | false
        'nf-core/'                  | false  // trailing slash
        '/nf-core/fastqc'           | false  // leading slash
    }
}
