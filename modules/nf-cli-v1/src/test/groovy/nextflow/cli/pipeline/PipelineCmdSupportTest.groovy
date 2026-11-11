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

import nextflow.config.RegistryConfig
import nextflow.exception.AbortOperationException
import spock.lang.Specification

class PipelineCmdSupportTest extends Specification {

    def 'should find the registry address of a name and hub' () {
        given:
        def config = new RegistryConfig(url: 'https://registry-dev.nextflow.io/api', apiKey: 'secret')

        when:
        def address = PipelineCmdSupport.registryAddress(NAME, HUB, config)

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
    }

    def 'should reject a hub that is neither a provider nor a registry' () {
        when:
        PipelineCmdSupport.registryAddress('nf-core/rnaseq', 'example.com', new RegistryConfig())

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Unknown repository provider or registry: `example.com`')
        e.message.contains('registries: registry.nextflow.io')
    }

    def 'should get the page URL of a registry' () {
        expect:
        PipelineCmdSupport.pageUrl('https://registry.nextflow.io/api') == 'https://registry.nextflow.io'
        PipelineCmdSupport.pageUrl('https://registry.nextflow.io/api/') == 'https://registry.nextflow.io'
        PipelineCmdSupport.pageUrl('https://example.com/registry') == 'https://example.com/registry'
    }
}
