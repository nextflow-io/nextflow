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

import java.nio.file.Files
import java.nio.file.Path

/**
 * Tests for ModuleManifest
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ModuleManifestTest extends Specification {

    @TempDir
    Path tempDir

    def 'should load valid manifest' () {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')
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
        def manifest = ModuleManifest.load(metaYaml)

        then:
        manifest.name == 'nf-core/fastqc'
        manifest.version == '1.0.0'
        manifest.description == 'FastQC quality control'
        manifest.authors == ['John Doe']
        manifest.license == 'MIT'
        manifest.keywords == ['quality-control', 'fastq']
        manifest.requires == ['nextflow': '>=24.04.0']
    }

    def 'should fail to load non-existent manifest' () {
        given:
        def metaYaml = tempDir.resolve('meta.yaml')

        when:
        ModuleManifest.load(metaYaml)

        then:
        thrown(AbortOperationException)
    }

    def 'should validate complete manifest' () {
        given:
        def manifest = new ModuleManifest(
            name: 'nf-core/fastqc',
            version: '1.0.0',
            description: 'FastQC quality control',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

        then:
        errors.isEmpty()
        manifest.isValid()
    }

    def 'should detect missing required fields' () {
        given:
        def manifest = new ModuleManifest(
            name: 'nf-core/fastqc'
            // missing version, description, license
        )

        when:
        def errors = manifest.validate()

        then:
        errors.size() == 3
        errors.any { it.contains('version') }
        errors.any { it.contains('description') }
        errors.any { it.contains('license') }
        !manifest.isValid()
    }

    def 'should validate version format' () {
        given:
        def manifest = new ModuleManifest(
            name: 'nf-core/fastqc',
            version: version,
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

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

    def 'should validate module name format' () {
        given:
        def manifest = new ModuleManifest(
            name: name,
            version: '1.0.0',
            description: 'Test',
            license: 'MIT'
        )

        when:
        def errors = manifest.validate()

        then:
        errors.isEmpty() == valid

        where:
        name                | valid
        'nf-core/fastqc'    | true
        'myorg/my-module'   | true
        'org_1/tool_2'      | true
        'fastqc'            | false
        '@nf-core/fastqc'   | false
        'nf-core/fast qc'   | false
    }
}
