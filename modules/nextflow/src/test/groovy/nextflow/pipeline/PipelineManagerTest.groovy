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

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

class PipelineManagerTest extends Specification {

    @TempDir
    Path base

    RegistryClient client = Mock(RegistryClient)
    PipelineGitFetcher fetcher = Mock(PipelineGitFetcher)

    PipelineManager manager() {
        new PipelineManager(base, client, fetcher)
    }

    static PipelineRelease release(String version, String commit = 'abc', PipelineReleaseStatus status = PipelineReleaseStatus.PUBLISHED) {
        new PipelineRelease().version(version).gitRevision(version).commitId(commit).status(status)
    }

    static Pipeline pipeline(PipelineRelease latest, boolean revoked = false) {
        new Pipeline().latest(latest).repositoryUrl('https://github.com/nf-core/rnaseq').workflowPath('main.nf').revoked(revoked)
    }

    def 'installs latest when missing'() {
        when:
        def main = manager().resolve('nf-core/rnaseq')
        then:
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline(release('3.14.0'))
        client.primaryUrl() >> 'https://reg.io'
        1 * fetcher.fetch('https://github.com/nf-core/rnaseq', '3.14.0', 'abc', _) >> { args -> args[3].resolve('main.nf').text = 'hello' }
        and:
        main == base.resolve('pipelines/nf-core/rnaseq/main.nf')
        def info = PipelineInfo.load(base.resolve('pipelines/nf-core/rnaseq'))
        info.version == '3.14.0'
        info.commitId == 'abc'
        info.source == 'git'
        info.registryUrl == 'https://reg.io'
        info.checksum
    }

    def 'null workflowPath defaults to main.nf'() {
        when:
        def info = manager().install('nf-core/rnaseq', null, false)
        then:
        1 * client.getPipeline(_) >> pipeline(release('1.0.0')).workflowPath(null)
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        info.workflowPath == 'main.nf'
    }

    def 'uses installed pipeline without calling the registry'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        m.resolve('nf-core/rnaseq')

        when:
        def main = m.resolve('nf-core/rnaseq')
        then:
        0 * client._
        0 * fetcher._
        main == base.resolve('pipelines/nf-core/rnaseq/main.nf')
    }

    def 'warns and keeps a modified pipeline'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        def main = m.resolve('nf-core/rnaseq')
        main.text = 'edited'

        when:
        def again = m.resolve('nf-core/rnaseq')
        then:
        0 * client._
        0 * fetcher._
        again == main
        main.text == 'edited'
    }

    def 'install refuses a modified pipeline without force, replaces it with force'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        def main = m.resolve('nf-core/rnaseq')
        main.text = 'edited'

        when:
        m.install('nf-core/rnaseq', null, false)
        then:
        thrown(AbortOperationException)
        main.text == 'edited'

        when:
        m.install('nf-core/rnaseq', null, true)
        then:
        main.text == 'x'
    }

    def 'install of the same version and checksum is a no-op'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        m.install('nf-core/rnaseq', null, false)

        when:
        def info = m.install('nf-core/rnaseq', null, false)
        then:
        0 * fetcher._
        info.version == '1.0.0'
    }

    def 'install -version replaces an unmodified different version'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('2.0.0'))
        client.getPipelineRelease('nf-core/rnaseq', '1.0.0') >> release('1.0.0', 'old')
        fetcher.fetch(_, _, 'old', _) >> { args -> args[3].resolve('main.nf').text = 'old' }
        fetcher.fetch(_, _, 'abc', _) >> { args -> args[3].resolve('main.nf').text = 'new' }
        m.install('nf-core/rnaseq', '1.0.0', false)

        when:
        def info = m.install('nf-core/rnaseq', null, false)
        then:
        info.version == '2.0.0'
        base.resolve('pipelines/nf-core/rnaseq/main.nf').text == 'new'
        // no leftover temp dirs
        Files.list(base.resolve('pipelines/nf-core')).count() == 1
    }

    def 'failed install leaves no directory'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))

        when:
        m.resolve('nf-core/rnaseq')
        then:
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('.git').toFile().mkdirs(); throw new AbortOperationException('boom') }
        thrown(AbortOperationException)
        !Files.exists(base.resolve('pipelines/nf-core/rnaseq'))
        Files.list(base.resolve('pipelines/nf-core')).count() == 0

        when:
        m.resolve('nf-core/rnaseq')
        then:
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        Files.exists(base.resolve('pipelines/nf-core/rnaseq/main.nf'))
    }

    def 'fails on no published release or revoked'() {
        when:
        manager().install('nf-core/rnaseq', null, false)
        then:
        client.getPipeline(_) >> new Pipeline().repositoryUrl('u')
        def e = thrown(AbortOperationException)
        e.message.contains('has no published release')

        when:
        manager().install('nf-core/rnaseq', null, false)
        then:
        client.getPipeline(_) >> pipeline(release('1.0.0'), true)
        e = thrown(AbortOperationException)
        e.message.contains('is revoked')

        when:
        manager().install('nf-core/rnaseq', '1.0.0', false)
        then:
        client.getPipeline(_) >> pipeline(release('2.0.0'))
        client.getPipelineRelease(_, _) >> release('1.0.0', 'abc', PipelineReleaseStatus.REVOKED)
        e = thrown(AbortOperationException)
        e.message.contains('is revoked')
    }

    def 'rejects invalid pipeline names'() {
        when:
        manager().resolve(name)
        then:
        0 * client._
        thrown(AbortOperationException)
        where:
        name << ['../x', 'nf-core/a/b', 'NF-core/x', 'x']
    }

    def 'list returns installed pipelines'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }

        expect:
        m.list() == []

        when:
        m.install('nf-core/rnaseq', null, false)
        m.install('acme/tool', null, false)
        then:
        m.list()*.name as Set == ['nf-core/rnaseq', 'acme/tool'] as Set
    }

    def 'remove deletes, refuses modified without force'() {
        given:
        def m = manager()
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> args[3].resolve('main.nf').text = 'x' }
        m.install('nf-core/rnaseq', null, false)
        def main = m.pipelineDir('nf-core/rnaseq').resolve('main.nf')
        main.text = 'edited'

        when:
        m.remove('nf-core/rnaseq', false)
        then:
        thrown(AbortOperationException)
        Files.exists(main)

        expect:
        m.remove('nf-core/rnaseq', true)
        !Files.exists(m.pipelineDir('nf-core/rnaseq'))
        !m.remove('nf-core/rnaseq', false)
    }
}
