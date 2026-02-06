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
import io.seqera.npr.api.schema.v1.ModuleSearchResult
import io.seqera.npr.api.schema.v1.SearchModulesResponse
import nextflow.exception.AbortOperationException
import nextflow.module.ModuleRegistryClient
import nextflow.cli.Launcher
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture


/**
 * Tests for ModuleSearch command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ModuleSearchTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    // No setup needed - using client field for mocking

    def 'should search and display results in formatted output'() {
        given:
        def result1 = new ModuleSearchResult(
            name: 'nf-core/fastqc',
            repositoryPath: 'nf-core/modules',
            description: 'FastQC quality control',
            relevanceScore: 0.95,
            keywords: ['quality-control', 'fastqc'],
            tools: ['fastqc'],
            revoked: false
        )
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
        def cmd = new ModuleSearch()
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
        def mockClient = Mock(ModuleRegistryClient) {
            search(_, _) >> response
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

    }

    def 'should search and display results in JSON output'() {
        given:
        def result1 = new ModuleSearchResult(
            name: 'nf-core/fastqc',
            description: 'FastQC quality control',
            relevanceScore: 0.95,
            keywords: ['quality-control'],
            tools: ['fastqc'],
            revoked: false
        )

        and:
        def cmd = new ModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['fastqc']
        cmd.limit = 10
        cmd.jsonOutput = true

        and:
        // Mock the registry client
        def mockClient = Mock(ModuleRegistryClient)
        mockClient.search('fastqc', 10) >> new SearchModulesResponse(
            query: 'fastqc',
            totalResults: 1,
            results: [result1]
        )
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString().readLines().last
        def json = new JsonSlurper().parseText(output)

        then:
        json.query == 'fastqc'
        json.totalResults == 1
        json.results.size() == 1
        json.results[0].name == 'nf-core/fastqc'
        json.results[0].description == 'FastQC quality control'

    }

    def 'should handle no search results'() {
        given:
        def cmd = new ModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['nonexistent-module']
        cmd.limit = 20

        and:
        // Mock empty results
        def mockClient = Mock(ModuleRegistryClient)
        mockClient.search('nonexistent-module', 20) >> new SearchModulesResponse(
            query: 'nonexistent-module',
            totalResults: 0,
            results: []
        )
        cmd.client = mockClient

        when:
        cmd.run()
        def output = capture.toString()

        then:
        output.contains('No modules found')

    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new ModuleSearch()
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
        def cmd = new ModuleSearch()
        cmd.launcher = Mock(Launcher){
            getOptions() >> null
        }
        cmd.args = ['test']
        cmd.limit = 5

        and:
        // Mock the client with 5 results
        def mockClient = Mock(ModuleRegistryClient)
        mockClient.search('test', 5) >> new SearchModulesResponse(
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
}
