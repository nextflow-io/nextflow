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
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.api.schema.v1.ModuleKind
import io.seqera.npr.api.schema.v1.ModuleSearchResult
import io.seqera.npr.api.schema.v1.ListPipelinesResponse
import io.seqera.npr.api.schema.v1.PipelineMetadata
import io.seqera.npr.api.schema.v1.PipelineSummary
import io.seqera.npr.api.schema.v1.SearchModulesResponse
import io.seqera.npr.client.RegistryException
import nextflow.exception.AbortOperationException
import nextflow.cli.Launcher
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture


/**
 * Tests for CmdModuleSearch command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CmdModuleSearchTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    // No setup needed - using client field for mocking

    def 'should search and display results in formatted output'() {
        given:
        def result1 = new ModuleSearchResult(
            name: 'nf-core/fastqc',
            repositoryPath: 'nf-core/modules',
            kind: ModuleKind.WORKFLOW,
            description: 'FastQC quality control',
            relevanceScore: 0.95,
            keywords: ['quality-control', 'fastqc'],
            tools: ['fastqc'],
            revoked: false
        )
        // result2 has no kind -> should default to Process
        def result2 = new ModuleSearchResult(
            name: 'nf-core/multiqc',
            repositoryPath: 'nf-core/modules',
            description: 'MultiQC reporting',
            relevanceScore: 0.85,
            keywords: ['quality-control', 'reporting'],
            tools: ['multiqc'],
            revoked: false
        )

        and:
        def cmd = new CmdModuleSearch()
        cmd.args = ['quality']
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.limit = 20

        and:
        def response = new SearchModulesResponse(
            query: 'quality',
            totalResults: 2,
            results: [result1, result2]
        )
        assert response.results
        // Mock the registry client directly using the test field
        def mockClient = Mock(RegistryClient) {
            searchModules(_, _) >> response
        }
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('Searching for')
        output.contains('quality')
        output.contains('nf-core/fastqc')
        output.contains('FastQC quality control')
        output.contains('nf-core/multiqc')
        output.contains('MultiQC reporting')

        and: 'the module kind is shown, defaulting to Process when absent'
        output.contains('Kind: Workflow')
        output.contains('Kind: Process')

    }

    def 'should search and display results in JSON output'() {
        given:
        def result1 = new ModuleSearchResult(
            name: 'nf-core/fastqc',
            kind: ModuleKind.WORKFLOW,
            description: 'FastQC quality control',
            relevanceScore: 0.95,
            keywords: ['quality-control'],
            tools: ['fastqc'],
            revoked: false
        )
        // result2 has no kind -> should default to Process in the JSON output
        def result2 = new ModuleSearchResult(
            name: 'nf-core/multiqc',
            description: 'MultiQC reporting',
            relevanceScore: 0.85,
            keywords: ['reporting'],
            tools: ['multiqc'],
            revoked: false
        )

        and:
        def cmd = new CmdModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['fastqc']
        cmd.limit = 10
        cmd.output = 'json'

        and:
        // Mock the registry client
        def mockClient = Mock(RegistryClient)
        mockClient.searchModules('fastqc', 10) >> new SearchModulesResponse(
            query: 'fastqc',
            totalResults: 2,
            results: [result1, result2]
        )
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString().readLines()
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('Searching for') ? line : null }.join("\n")
        def json = new JsonSlurper().parseText(output)

        then:
        json.query == 'fastqc'
        json.totalResults == 2
        json.results.size() == 2
        json.results[0].name == 'nf-core/fastqc'
        json.results[0].description == 'FastQC quality control'

        and: 'kind is emitted, defaulting to Process when absent'
        json.results[0].kind == 'Workflow'
        json.results[1].kind == 'Process'

    }

    def 'should handle no search results'() {
        given:
        def cmd = new CmdModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['nonexistent-module']
        cmd.limit = 20

        and:
        // Mock empty results
        def mockClient = Mock(RegistryClient)
        mockClient.searchModules('nonexistent-module', 20) >> new SearchModulesResponse(
            query: 'nonexistent-module',
            totalResults: 0,
            results: []
        )
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('No modules or pipelines found')

    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new CmdModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = []

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should handle search with custom limit'() {
        given:
        def results = (1..5).collect { i ->
            new ModuleSearchResult(
                name: "nf-core/module${i}",
                description: "Module ${i}",
                relevanceScore: 0.9 - (i * 0.1),
                keywords: ['test'],
                tools: ["tool${i}"],
                revoked: false
            )
        }

        and:
        def cmd = new CmdModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['test']
        cmd.limit = 5

        and:
        // Mock the client with 5 results
        def mockClient = Mock(RegistryClient)
        mockClient.searchModules('test', 5) >> new SearchModulesResponse(
            query: 'test',
            totalResults: 5,
            results: results
        )
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('nf-core/module1')
        output.contains('nf-core/module5')
        (1..5).every { i -> output.contains("Module ${i}") }
    }

    private CmdModuleSearch searchCmd(RegistryClient client, String output = 'simple') {
        def cmd = new CmdModuleSearch(args: ['rna'], limit: 5, output: output, client: client)
        cmd.launcher = Mock(Launcher) { getOptions() >> null }
        return cmd
    }

    private static SearchModulesResponse modules(String... names) {
        new SearchModulesResponse(query: 'rna', totalResults: names.size(), results: names.collect { new ModuleSearchResult(name: it, description: "Module $it") })
    }

    private static ListPipelinesResponse pipelines(String... names) {
        new ListPipelinesResponse().pipelines(names.collect { new PipelineSummary().name(it).latestVersion('3.14.0').metadata(new PipelineMetadata().description('RNA-seq analysis')) })
    }

    def 'should list modules first, then pipelines'() {
        given:
        def client = Mock(RegistryClient)

        when:
        searchCmd(client).run()
        def output = capture.toString()

        then:
        1 * client.searchModules('rna', 5) >> modules('nf-core/salmon')
        1 * client.listPipelines('rna', 1, 4) >> pipelines('nf-core/rnaseq')
        output.contains('Top 2 matching result(s):')
        output.indexOf('nf-core/salmon') < output.indexOf('nf-core/rnaseq')
        output.contains('  nf-core/rnaseq\n    Kind: Pipeline\n    Version: 3.14.0\n    Description: RNA-seq analysis')
    }

    def 'should append pipelines to the JSON results'() {
        given:
        def client = Mock(RegistryClient) {
            searchModules('rna', 5) >> new SearchModulesResponse(query: 'rna', totalResults: 40, results: [new ModuleSearchResult(name: 'nf-core/salmon')])
            listPipelines('rna', 1, 4) >> pipelines('nf-core/rnaseq').totalSize(7L)
        }

        when:
        searchCmd(client, 'json').run()
        def output = capture.toString()
        def json = new JsonSlurper().parseText(output.substring(output.indexOf('{')))

        then: 'the registry totals of both kinds'
        json.totalResults == 47
        json.results*.name == ['nf-core/salmon', 'nf-core/rnaseq']
        json.results[0].kind == 'Process'
        json.results[1].kind == 'Pipeline'
        json.results[1].latestVersion == '3.14.0'
        json.results[1].description == 'RNA-seq analysis'
    }

    def 'should show pipelines when no module matches'() {
        given:
        def client = Mock(RegistryClient) {
            searchModules('rna', 5) >> modules()
            listPipelines('rna', 1, 5) >> pipelines('nf-core/rnaseq')
        }

        when:
        searchCmd(client).run()
        def output = capture.toString()

        then:
        output.contains('nf-core/rnaseq')
        !output.contains('No modules or pipelines found')
    }

    def 'should not search pipelines when modules fill the limit'() {
        given:
        def client = Mock(RegistryClient)

        when:
        searchCmd(client, 'json').run()
        def output = capture.toString()
        def json = new JsonSlurper().parseText(output.substring(output.indexOf('{')))

        then:
        1 * client.searchModules('rna', 5) >> new SearchModulesResponse(query: 'rna', totalResults: 12, results: (1..5).collect { new ModuleSearchResult(name: "nf-core/m$it") })
        0 * client.listPipelines(*_)
        json.results.size() == 5
        json.totalResults == 12
    }

    def 'should read further pipeline pages up to the limit'() {
        given:
        def client = Mock(RegistryClient)
        def cmd = searchCmd(client)
        cmd.limit = 150
        def page = { int n -> pipelines(*(1..n).collect { "nf-core/p$it" } as String[]) }

        when:
        cmd.run()

        then:
        1 * client.searchModules('rna', 150) >> modules()
        1 * client.listPipelines('rna', 1, 100) >> page(100).last(false)
        1 * client.listPipelines('rna', 2, 100) >> page(100).last(false)
        0 * client.listPipelines(*_)
        capture.toString().contains('Top 150 matching result(s):')
    }

    def 'should show modules only when the pipeline search fails'() {
        given:
        def client = Mock(RegistryClient) {
            searchModules('rna', 5) >> modules('nf-core/salmon')
            listPipelines('rna', 1, 4) >> { throw new RegistryException('Server error', 500) }
        }

        when:
        searchCmd(client).run()
        def output = capture.toString()

        then:
        output.contains('Top 1 matching result(s):')
        output.contains('nf-core/salmon')
    }

    def 'should report when neither modules nor pipelines match'() {
        given:
        def client = Mock(RegistryClient) {
            searchModules('rna', 5) >> modules()
            listPipelines('rna', 1, 5) >> new ListPipelinesResponse().pipelines([])
        }

        when:
        searchCmd(client).run()

        then:
        capture.toString().contains('No modules or pipelines found')
    }
}
