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

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.exception.AbortOperationException
import spock.lang.Specification

class PipelineSourceTest extends Specification {

    RegistryClient client = Mock(RegistryClient)

    def 'should resolve a registry release'() {
        given:
        def latest = new PipelineRelease().version('3.14.0').gitRevision('3.14.0').commitId('abc123')
        def pipeline = new Pipeline().repositoryUrl('https://github.com/nf-core/rnaseq').workflowPath('main.nf').latest(latest)
        client.primaryUrl() >> 'https://reg.io'

        when:
        def result = PipelineSource.resolve(client, 'nf-core/rnaseq', null)

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline
        0 * client.getPipelineRelease(_, _)
        result == new PipelineSource(
            name: 'nf-core/rnaseq',
            version: '3.14.0',
            workflowPath: 'main.nf',
            repositoryUrl: 'https://github.com/nf-core/rnaseq',
            revision: '3.14.0',
            commitId: 'abc123',
            registryUrl: 'https://reg.io')

        when:
        result = PipelineSource.resolve(client, 'nf-core/rnaseq', '3.13.0')

        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline
        1 * client.getPipelineRelease('nf-core/rnaseq', '3.13.0') >> new PipelineRelease().version('3.13.0').commitId('def456')
        result.version == '3.13.0'
        result.commitId == 'def456'
    }

    def 'should fail when the registry cannot resolve the release'() {
        given:
        def pipeline = new Pipeline().repositoryUrl('https://github.com/nf-core/rnaseq')
            .latest(LATEST ? new PipelineRelease().version('3.14.0').commitId('abc123').status(STATUS) : null)
            .revoked(REVOKED)
        client.getPipeline('nf-core/rnaseq') >> { if( PIPELINE_ERROR ) throw PIPELINE_ERROR; pipeline }
        client.getPipelineRelease('nf-core/rnaseq', 'dev') >> { throw new RegistryException('Not found', 404) }

        when:
        PipelineSource.resolve(client, NAME, VERSION)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains(MESSAGE)

        where:
        NAME             | VERSION | PIPELINE_ERROR                          | LATEST | STATUS                        | REVOKED | MESSAGE
        'hello'          | null    | null                                    | true   | null                          | false   | 'Invalid pipeline name'
        'nf-core/a/b'    | null    | null                                    | true   | null                          | false   | 'Invalid pipeline name'
        'nf-core/rnaseq' | null    | new RegistryException('Not found', 404) | true   | null                          | false   | '`nf-core/rnaseq` not found in the registry -- Not found'
        'nf-core/rnaseq' | 'dev'   | null                                    | true   | null                          | false   | '`nf-core/rnaseq@dev` not found in the registry -- Not found'
        'nf-core/rnaseq' | null    | null                                    | false  | null                          | false   | 'has no published release'
        'nf-core/rnaseq' | null    | null                                    | true   | PipelineReleaseStatus.REVOKED | false   | 'is revoked'
        'nf-core/rnaseq' | null    | null                                    | true   | null                          | true    | 'is revoked'
    }

    def 'should rethrow registry errors other than not found'() {
        given:
        client.getPipeline(_) >> { throw new RegistryException('Unavailable', 503) }

        when:
        PipelineSource.resolve(client, 'nf-core/rnaseq', null)

        then:
        def e = thrown(RegistryException)
        e.statusCode == 503
    }
}
