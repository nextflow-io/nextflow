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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.AbortOperationException
import nextflow.pipeline.PipelineInfo
import spock.lang.Specification
import spock.lang.TempDir

class CmdPipelineRemoveTest extends Specification {

    @TempDir
    Path base

    Path seed() {
        def dir = base.resolve('pipelines/nf-core/rnaseq')
        Files.createDirectories(dir)
        Files.writeString(dir.resolve('main.nf'), 'x')
        new PipelineInfo(name: 'nf-core/rnaseq', version: '3.14.0', checksum: PipelineInfo.checksum(dir)).save(dir)
        return dir
    }

    def 'should remove a pipeline'() {
        given:
        def dir = seed()

        when:
        new CmdPipelineRemove(args: ['nf-core/rnaseq'], root: base).run()

        then:
        !Files.exists(dir)
    }

    def 'should refuse to remove a modified copy without -force'() {
        given:
        def dir = seed()
        Files.writeString(dir.resolve('main.nf'), 'changed')

        when:
        new CmdPipelineRemove(args: ['nf-core/rnaseq'], root: base).run()

        then:
        thrown(AbortOperationException)
        Files.exists(dir)

        when:
        new CmdPipelineRemove(args: ['nf-core/rnaseq'], root: base, force: true).run()

        then:
        !Files.exists(dir)
    }
}
