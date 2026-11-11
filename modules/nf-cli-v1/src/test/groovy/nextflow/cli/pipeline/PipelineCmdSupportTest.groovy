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

import java.nio.file.Path

import nextflow.cli.CmdBase
import nextflow.cli.Launcher
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

class PipelineCmdSupportTest extends Specification {

    @TempDir
    Path folder

    CmdBase cmd() {
        def cmd = Mock(CmdBase)
        cmd.getLauncher() >> Mock(Launcher) { getOptions() >> null }
        return cmd
    }

    def 'should find the registry address of a name and hub' () {
        given:
        def config = folder.resolve('registry.config')
        config.text = "registry { url = 'https://registry-dev.nextflow.io/api'; apiKey = 'secret' }"

        when:
        def address = PipelineCmdSupport.registryAddress(cmd(), NAME, HUB, [config])

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
    }

    def 'should reject a hub that is neither a provider nor a registry' () {
        when:
        PipelineCmdSupport.registryAddress(cmd(), 'nf-core/rnaseq', 'example.com')

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Unknown repository provider or registry: `example.com`')
        e.message.contains('registries: registry.nextflow.io')
    }
}
