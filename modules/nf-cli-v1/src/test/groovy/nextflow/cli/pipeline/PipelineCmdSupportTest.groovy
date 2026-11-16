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

package nextflow.cli.pipeline

import com.beust.jcommander.ParameterException
import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import spock.lang.Specification

class PipelineCmdSupportTest extends Specification {

    def 'should find the registry address of a name and hub' () {
        given:
        def config = new RegistryConfig(url: 'https://registry-dev.nextflow.io/api', apiKey: 'secret')

        when:
        def address = PipelineCmdSupport.registryAddress(NAME, HUB) { config }

        then:
        address?.url == API
        address?.name == PIPELINE
        address?.apiKey == KEY

        where:
        NAME                                                         | HUB                        | API                                     | PIPELINE         | KEY
        'nf-core/rnaseq'                                             | null                       | null                                    | null             | null
        'nf-core/rnaseq'                                             | 'github'                   | null                                    | null             | null
        'https://github.com/nf-core/rnaseq'                          | null                       | null                                    | null             | null
        'nf-core/rnaseq'                                             | 'registry-dev.nextflow.io' | 'https://registry-dev.nextflow.io/api'  | 'nf-core/rnaseq' | 'secret'
        'https://registry-dev.nextflow.io/pipelines/nf-core/rnaseq'  | null                       | 'https://registry-dev.nextflow.io/api'  | 'nf-core/rnaseq' | 'secret'
        'https://registry.nextflow.io/pipelines/nf-core/rnaseq'      | null                       | 'https://registry.nextflow.io/api'      | 'nf-core/rnaseq' | null
        'nf-core/rnaseq'                                             | 'registry.nextflow.io'     | 'https://registry.nextflow.io/api'      | 'nf-core/rnaseq' | null
        'nf-core/rnaseq'                                             | 'Registry-Dev.Nextflow.io' | 'https://registry-dev.nextflow.io/api'  | 'nf-core/rnaseq' | 'secret'
        'https://registry-dev.nextflow.io/pipelines/acme/hello/'     | null                       | 'https://registry-dev.nextflow.io/api'  | 'acme/hello'     | 'secret'
        'https://registry.nextflow.io/modules/nf-core/fastqc'        | null                       | null                                    | null             | null
        'https://example.com/pipelines/nf-core/rnaseq'               | null                       | null                                    | null             | null
        'http://registry-dev.nextflow.io/pipelines/nf-core/rnaseq'   | null                       | 'https://registry-dev.nextflow.io/api'  | 'nf-core/rnaseq' | 'secret'
    }

    def 'should reject a version in a registry address' () {
        when:
        PipelineCmdSupport.registryAddress(NAME, HUB) { new RegistryConfig() }

        then:
        def e = thrown(AbortOperationException)
        e.message == "Not a valid registry pipeline address: nf-core/rnaseq@3.14.0 -- give the version with the `-r` option"

        where:
        NAME                                                            | HUB
        'https://registry.nextflow.io/pipelines/nf-core/rnaseq@3.14.0'  | null
        'nf-core/rnaseq@3.14.0'                                         | 'registry.nextflow.io'
    }

    def 'should match a registry URL under a base path' () {
        given:
        def config = new RegistryConfig(url: 'https://example.com/registry/api')

        expect:
        PipelineCmdSupport.registryAddress('https://example.com/registry/pipelines/acme/hello', null) { config }?.name == 'acme/hello'
        PipelineCmdSupport.registryAddress('https://example.com/pipelines/acme/hello', null) { config } == null
    }

    def 'should reject a hub that is neither a provider nor a registry' () {
        when:
        PipelineCmdSupport.registryAddress('nf-core/rnaseq', 'example.com') { new RegistryConfig() }

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Unknown repository provider or registry: `example.com`')
        e.message.contains('registries: registry.nextflow.io')
    }

    def 'should not read the config for a Git address' () {
        expect:
        PipelineCmdSupport.registryAddress(NAME, HUB) { throw new AssertionError('config must not be read') } == null

        where:
        NAME                                       | HUB
        'nf-core/rnaseq'                           | null
        'nf-core/rnaseq'                           | 'github'
        'nf-core/rnaseq'                           | 'gitlab'
        'https://github.com/nf-core/rnaseq'        | null
        'https://github.com/nf-core/rnaseq'        | 'registry.nextflow.io'
        'https://host/pipelines/a b'               | null
        'https://gitlab.com/acme/pipelines/rnaseq' | null
        'https://GitHub.com/acme/pipelines/rnaseq' | null
        'https://bitbucket.org/acme/pipelines'     | null
    }

    def 'should match a registry by host and port' () {
        given:
        def config = new RegistryConfig(url: 'http://localhost:8080/api')

        expect:
        PipelineCmdSupport.registryAddress(NAME, HUB) { config }?.url == 'http://localhost:8080/api'

        where:
        NAME                                                    | HUB
        'acme/hello'                                            | 'localhost'
        'acme/hello'                                            | 'localhost:8080'
        'http://localhost:8080/pipelines/acme/hello'            | null
    }

    def 'should get the page URL of a registry' () {
        expect:
        PipelineCmdSupport.pageUrl('https://registry.nextflow.io/api') == 'https://registry.nextflow.io'
        PipelineCmdSupport.pageUrl('https://registry.nextflow.io/api/') == 'https://registry.nextflow.io'
        PipelineCmdSupport.pageUrl('https://example.com/registry') == 'https://example.com/registry'
    }

    def 'should accept the output modes of each command' () {
        when:
        VALIDATOR.getDeclaredConstructor().newInstance().validate('-o', MODE)

        then:
        noExceptionThrown()

        where:
        VALIDATOR                                  | MODE
        PipelineCmdSupport.TextOutputValidator     | 'text'
        PipelineCmdSupport.TextOutputValidator     | 'json'
        PipelineCmdSupport.SimpleOutputValidator   | 'simple'
        PipelineCmdSupport.SimpleOutputValidator   | 'json'
    }

    def 'should reject the output modes of other commands' () {
        when:
        VALIDATOR.getDeclaredConstructor().newInstance().validate('-o', MODE)

        then:
        def e = thrown(ParameterException)
        e.message == "Output mode must be one of $EXPECTED (found: $MODE)".toString()

        where:
        VALIDATOR                                  | MODE     | EXPECTED
        PipelineCmdSupport.TextOutputValidator     | 'simple' | ['text', 'json']
        PipelineCmdSupport.SimpleOutputValidator   | 'text'   | ['simple', 'json']
    }
}
