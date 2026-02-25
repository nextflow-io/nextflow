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

package nextflow.cli.module

import groovy.json.JsonSlurper
import io.seqera.npr.api.schema.v1.Module
import io.seqera.npr.api.schema.v1.ModuleChannel
import io.seqera.npr.api.schema.v1.ModuleChannelItem
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.api.schema.v1.ModuleRelease
import io.seqera.npr.api.schema.v1.ModuleTool
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleRegistryClient
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture

import java.nio.file.Path

/**
 * Tests for ModuleInfo command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleInfoTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    @TempDir
    Path tempDir

    def 'should display module info in formatted output'() {
        given:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control analysis',
            authors: ['nf-core', 'community'],
            keywords: ['quality-control', 'fastqc', 'reads']
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            description: 'FastQC module',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Module:')
        output.contains('nf-core/fastqc')
        output.contains('Version:')
        output.contains('1.0.0')
        output.contains('Description:')
        output.contains('FastQC quality control analysis')
        output.contains('Authors:')
        output.contains('nf-core, community')
        output.contains('Keywords:')
        output.contains('quality-control, fastqc, reads')
        output.contains('Usage Template:')
    }

    def 'should display module info with specific version'() {
        given:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control'
        )

        and:
        def release = new ModuleRelease(
            version: '0.9.0',
            description: 'FastQC module',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.version = '0.9.0'
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockClient = Mock(ModuleRegistryClient) {
            fetchRelease(_, _) >> release
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Version:')
        output.contains('0.9.0')
    }

    def 'should display module info in JSON format'() {
        given:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            authors: ['nf-core'],
            keywords: ['quality-control']
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            description: 'FastQC module',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.jsonOutput = true
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()
        // Extract JSON part (skip debug/log lines)
        def lines = output.readLines()
        def jsonStart = lines.findIndexOf { it.trim().startsWith('{') }
        def jsonText = lines[jsonStart..-1].join('\n')
        def json = new JsonSlurper().parseText(jsonText)

        then:
        json.name == 'nf-core/fastqc'
        json.fullName == '@nf-core/fastqc'
        json.version == '1.0.0'
        json.description == 'FastQC quality control'
        json.authors == ['nf-core']
        json.keywords == ['quality-control']
        json.usageTemplate != null
    }

    def 'should display module info with tools'() {
        given:
        def tool = new ModuleTool(
            name: 'fastqc',
            version: '0.12.1',
            homepage: URI.create('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'),
            documentation: URI.create('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/')
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            tools: [tool]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Tools:')
        output.contains('fastqc v0.12.1')
        output.contains('Homepage: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/')
    }

    def 'should display module info with inputs'() {
        given:
        def inputItem1 = new ModuleChannelItem(
            name: 'reads',
            type: 'file',
            description: 'Input FASTQ files',
            pattern: '*.fastq.gz'
        )
        def inputItem2 = new ModuleChannelItem(
            name: 'meta',
            type: 'map',
            description: 'Sample metadata'
        )

        and:
        def input = new ModuleChannel(
            tuple: true,
            items: [inputItem1, inputItem2]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            input: [input]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Input:')
        output.contains('(tuple)')
        output.contains('reads (file)')
        output.contains('Input FASTQ files')
        output.contains('Pattern: *.fastq.gz')
        output.contains('meta (map)')
        output.contains('Sample metadata')
    }

    def 'should display module info with outputs'() {
        given:
        def outputItem = new ModuleChannelItem(
            name: 'html',
            type: 'file',
            description: 'FastQC HTML report',
            pattern: '*_fastqc.html'
        )

        and:
        def outputChannel = new ModuleChannel(
            tuple: false,
            items: [outputItem]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            output: ['html': outputChannel]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Output:')
        output.contains('html')
        output.contains('html (file)')
        output.contains('FastQC HTML report')
        output.contains('Pattern: *_fastqc.html')
    }

    def 'should generate usage template with inputs'() {
        given:
        def inputItem1 = new ModuleChannelItem(
            name: 'reads',
            type: 'file'
        )
        def inputItem2 = new ModuleChannelItem(
            name: 'sample-id',
            type: 'val'
        )

        and:
        def input1 = new ModuleChannel(
            tuple: false,
            items: [inputItem1]
        )
        def input2 = new ModuleChannel(
            tuple: false,
            items: [inputItem2]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            input: [input1, input2]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Usage Template:')
        output.contains('nextflow module run nf-core/fastqc')
        output.contains('--reads <FASTQ_FILE>')
        output.contains('--sample-id <SAMPLE_ID>')
    }

    def 'should generate usage template with version'() {
        given:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            input: []
        )

        and:
        def release = new ModuleRelease(
            version: '2.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.version = '2.0.0'
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockClient = Mock(ModuleRegistryClient) {
            fetchRelease(_, _) >> release
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Usage Template:')
        output.contains('nextflow module run nf-core/fastqc')
        output.contains('-version 2.0.0')
    }

    def 'should generate usage template with map inputs'() {
        given:
        def inputItem = new ModuleChannelItem(
            name: 'params',
            type: 'map'
        )

        and:
        def input = new ModuleChannel(
            tuple: false,
            items: [inputItem]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            input: [input]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Usage Template:')
        output.contains('--params.<KEY> <PARAMS_KEY_VALUE>')
    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new ModuleInfo()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = []
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail with multiple arguments'() {
        given:
        def cmd = new ModuleInfo()
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.args = ['module1', 'module2']
        cmd.root = tempDir

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should display minimal info when metadata is sparse'() {
        given:
        def metadata = new ModuleMetadata(
            description: null,
            authors: null,
            keywords: null
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            description: 'Module description',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Module:')
        output.contains('nf-core/fastqc')
        output.contains('Version:')
        output.contains('1.0.0')
        output.contains('Description:')
        output.contains('Module description')
    }

    def 'should display N/A when no description is available'() {
        given:
        def metadata = new ModuleMetadata(
            description: null
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            description: null,
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Description: N/A')
    }

    def 'should display module with tuple outputs'() {
        given:
        def outputItem1 = new ModuleChannelItem(
            name: 'html',
            type: 'file',
            description: 'HTML report'
        )
        def outputItem2 = new ModuleChannelItem(
            name: 'zip',
            type: 'file',
            description: 'ZIP archive'
        )

        and:
        def outputChannel = new ModuleChannel(
            tuple: true,
            items: [outputItem1, outputItem2]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            output: ['qc': outputChannel]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Output:')
        output.contains('qc (tuple)')
        output.contains('html (file)')
        output.contains('HTML report')
        output.contains('zip (file)')
        output.contains('ZIP archive')
        output.contains(')')
    }

    def 'should include all tool information in JSON output'() {
        given:
        def tool = new ModuleTool(
            name: 'fastqc',
            version: '0.12.1',
            homepage: URI.create('https://example.com'),
            documentation: URI.create('https://docs.example.com')
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            tools: [tool]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.jsonOutput = true
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()
        // Extract JSON part (skip debug/log lines)
        def lines = output.readLines()
        def jsonStart = lines.findIndexOf { it.trim().startsWith('{') }
        def jsonText = lines[jsonStart..-1].join('\n')
        def json = new JsonSlurper().parseText(jsonText)

        then:
        json.tools.size() == 1
        json.tools[0].name == 'fastqc'
        json.tools[0].version == '0.12.1'
        json.tools[0].homepage.scheme == 'https'
        json.tools[0].homepage.host == 'example.com'
        json.tools[0].documentation.scheme == 'https'
        json.tools[0].documentation.host == 'docs.example.com'
    }

    def 'should include input/output information in JSON output'() {
        given:
        def inputItem = new ModuleChannelItem(
            name: 'reads',
            type: 'file',
            description: 'Input reads',
            pattern: '*.fastq.gz'
        )
        def inputChannel = new ModuleChannel(
            tuple: true,
            items: [inputItem]
        )

        and:
        def outputItem = new ModuleChannelItem(
            name: 'html',
            type: 'file',
            description: 'HTML report',
            pattern: '*.html'
        )
        def outputChannel = new ModuleChannel(
            tuple: false,
            items: [outputItem]
        )

        and:
        def metadata = new ModuleMetadata(
            description: 'FastQC quality control',
            input: [inputChannel],
            output: ['html': outputChannel]
        )

        and:
        def release = new ModuleRelease(
            version: '1.0.0',
            metadata: metadata
        )

        and:
        def cmd = new ModuleInfo()
        cmd.args = ['nf-core/fastqc']
        cmd.jsonOutput = true
        cmd.launcher = Mock(Launcher) {
            getOptions() >> null
        }
        cmd.root = tempDir

        and:
        def mockModule = Stub(Module) {
            getLatest() >> release
        }
        def mockClient = Mock(ModuleRegistryClient) {
            fetchModule(_) >> mockModule
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()
        // Extract JSON part (skip debug/log lines)
        def lines = output.readLines()
        def jsonStart = lines.findIndexOf { it.trim().startsWith('{') }
        def jsonText = lines[jsonStart..-1].join('\n')
        def json = new JsonSlurper().parseText(jsonText)

        then:
        json.input.size() == 1
        json.input[0].tuple == true
        json.input[0].items[0].name == 'reads'
        json.input[0].items[0].type == 'file'
        json.input[0].items[0].description == 'Input reads'
        json.input[0].items[0].pattern == '*.fastq.gz'

        and:
        json.output.html.tuple == false
        json.output.html.items[0].name == 'html'
        json.output.html.items[0].type == 'file'
        json.output.html.items[0].description == 'HTML report'
        json.output.html.items[0].pattern == '*.html'
    }
}