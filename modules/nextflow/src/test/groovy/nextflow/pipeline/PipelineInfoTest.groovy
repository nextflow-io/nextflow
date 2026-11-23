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

package nextflow.pipeline

import java.nio.file.Files
import java.nio.file.Path

import spock.lang.Specification
import spock.lang.TempDir

class PipelineInfoTest extends Specification {

    @TempDir
    Path dir

    def 'save then load round-trips every field'() {
        given:
        def info = new PipelineInfo(name: 'nf-core/rnaseq', version: '1.0.0', workflowPath: 'main.nf', source: 'artifact',
                repositoryUrl: 'https://github.com/nf-core/rnaseq', revision: 'v1', commitId: 'abc', checksum: 'dead', registryUrl: 'https://r.io',
                artifactUrl: 'https://a.io/b.tgz', artifactChecksum: 'sha256:beef')
        when:
        info.save(dir)
        def loaded = PipelineInfo.load(dir)
        then:
        loaded.properties.findAll { it.key != 'class' } == info.properties.findAll { it.key != 'class' }
        loaded.registryUrl == 'https://r.io'
        loaded.artifactUrl == 'https://a.io/b.tgz'
        loaded.artifactChecksum == 'sha256:beef'
    }

    def 'load returns null when .pipeline-info is missing'() {
        expect:
        PipelineInfo.load(dir) == null
    }

    def 'checksum is stable and changes when a file changes'() {
        given:
        def main = dir.resolve('main.nf')
        main.text = 'one'
        def c1 = PipelineInfo.checksum(dir)
        expect:
        c1 ==~ /[0-9a-f]{64}/
        PipelineInfo.checksum(dir) == c1
        when:
        main.text = 'two'
        then:
        PipelineInfo.checksum(dir) != c1
    }

    def 'checksum ignores .pipeline-info and .git'() {
        given:
        dir.resolve('main.nf').text = 'one'
        def c1 = PipelineInfo.checksum(dir)
        when:
        new PipelineInfo(name: 'x').save(dir)
        Files.createDirectories(dir.resolve('.git'))
        dir.resolve('.git/HEAD').text = 'ref'
        then:
        PipelineInfo.checksum(dir) == c1
        and:
        !new PipelineInfo(checksum: c1).isModified(dir)
    }

    def 'ignores nested dirs with their own info file'() {
        given:
        dir.resolve('main.nf').text = 'one'
        def c1 = PipelineInfo.checksum(dir)
        when:
        def mod = Files.createDirectories(dir.resolve('modules/nf-core/fastqc'))
        mod.resolve('main.nf').text = 'm'
        mod.resolve('.module-info').text = 'x=1'
        def pipe = Files.createDirectories(dir.resolve('pipelines/x/y'))
        pipe.resolve('main.nf').text = 'p'
        pipe.resolve('.pipeline-info').text = 'x=1'
        then:
        PipelineInfo.checksum(dir) == c1
        when:
        def other = Files.createDirectories(dir.resolve('modules/nf-core/other'))
        other.resolve('main.nf').text = 'o'
        then:
        PipelineInfo.checksum(dir) != c1
    }
}
