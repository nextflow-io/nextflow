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

import java.nio.file.Path

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.config.RegistryConfig
import nextflow.module.spi.RemoteModuleResolverProvider
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Test for DefaultRemoteModuleResolver SPI implementation
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultRemoteModuleResolverTest extends Specification {

    def 'should load resolver via SPI'() {
        when:
        def resolver = RemoteModuleResolverProvider.getInstance()

        then:
        resolver != null
        resolver.class.name == 'nextflow.module.DefaultRemoteModuleResolver'
        resolver.priority == 0
    }

    def 'should return default priority'() {
        given:
        def resolver = new DefaultRemoteModuleResolver()

        expect:
        resolver.getPriority() == 0
    }

    @TempDir
    Path base

    def 'should auto-install a pipeline'() {
        given:
        def client = Mock(RegistryClient)
        def fetcher = Mock(PipelineGitFetcher)
        def resolver = new DefaultRemoteModuleResolver() {
            @Override
            protected ModuleResolver createResolver(Path baseDir, RegistryConfig registryConfig) {
                return new ModuleResolver(baseDir, client, fetcher)
            }
        }

        when:
        def main = resolver.resolve('nf-core/rnaseq', base)

        then:
        1 * client.getModule('nf-core/rnaseq') >> { throw new RegistryException('Not found', 404) }
        1 * client.getPipeline('nf-core/rnaseq') >> new Pipeline()
            .repositoryUrl('https://github.com/nf-core/rnaseq')
            .latest(new PipelineRelease().version('3.14.0').commitId('abc'))
        1 * fetcher.fetch('https://github.com/nf-core/rnaseq', null, 'abc', _) >> { args -> args[3].resolve('main.nf').text = 'workflow { }'; null }
        and:
        main == base.resolve('modules/nf-core/rnaseq/main.nf')
        ModuleInfo.load(main.parent, 'kind') == ModuleInfo.KIND_PIPELINE
    }
}
