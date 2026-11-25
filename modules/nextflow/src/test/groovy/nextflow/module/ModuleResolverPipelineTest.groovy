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

import java.nio.file.Files
import java.nio.file.Path

import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.exception.AbortOperationException
import spock.lang.Specification
import spock.lang.TempDir

/**
 * Tests for installing registry pipelines with {@link ModuleResolver}
 */
class ModuleResolverPipelineTest extends Specification {

    @TempDir
    Path base

    RegistryClient client = Mock(RegistryClient)
    PipelineGitFetcher fetcher = Mock(PipelineGitFetcher)

    def ref = new ModuleReference('nf-core', 'rnaseq')

    ModuleResolver resolver() {
        new ModuleResolver(base, client, fetcher)
    }

    Path dir() {
        base.resolve('modules/nf-core/rnaseq')
    }

    static RegistryException notFound() {
        new RegistryException('Not found', 404)
    }

    static PipelineRelease release(String version, String commit = 'abc', PipelineReleaseStatus status = PipelineReleaseStatus.PUBLISHED) {
        new PipelineRelease().version(version).gitRevision(version).commitId(commit).status(status)
    }

    static Pipeline pipeline(PipelineRelease latest, boolean revoked = false) {
        new Pipeline().latest(latest).repositoryUrl('https://github.com/nf-core/rnaseq').workflowPath('main.nf').revoked(revoked)
    }

    static String writeMain(Path target, String text, String commit = null) {
        target.resolve('main.nf').text = text
        return commit
    }

    def 'installs a pipeline when missing'() {
        when:
        def main = resolver().resolve(ref, null, true)
        then:
        1 * client.getModule('nf-core/rnaseq') >> { throw notFound() }
        1 * client.getPipeline('nf-core/rnaseq') >> pipeline(release('3.14.0'))
        client.primaryUrl() >> 'https://reg.io'
        1 * fetcher.fetch('https://github.com/nf-core/rnaseq', '3.14.0', 'abc', _) >> { args -> writeMain(args[3] as Path, 'hello') }
        and:
        main == dir().resolve('main.nf')
        def info = ModuleInfo.load(dir())
        info.kind == 'Pipeline'
        info.version == '3.14.0'
        info.repositoryUrl == 'https://github.com/nf-core/rnaseq'
        info.revision == '3.14.0'
        info.commitId == 'abc'
        info.workflowPath == 'main.nf'
        info.registryUrl == 'https://reg.io'
        info.checksum
        and:
        def installed = new ModuleStorage(base).getInstalledModule(ref)
        installed.pipeline
        installed.installedVersion == '3.14.0'
        installed.integrity == ModuleIntegrity.VALID
        Files.list(dir().parent).count() == 1
    }

    def 'records the fetched commit when the release has no commit id'() {
        when:
        resolver().installWithDependencies(ref)
        then:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0', null))
        1 * fetcher.fetch(_, '1.0.0', null, _) >> { args -> writeMain(args[3] as Path, 'x', 'def') }
        ModuleInfo.load(dir(), 'commitId') == 'def'
    }

    def 'reuses an installed pipeline without a registry call'() {
        given:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> writeMain(args[3] as Path, 'x') }
        resolver().installWithDependencies(ref)

        when:
        def main = resolver().resolve(ref, null, true)
        resolver().installWithDependencies(ref, '1.0.0')
        then:
        0 * client._
        0 * fetcher._
        main == dir().resolve('main.nf')
    }

    def 'refuses to replace a modified pipeline, replaces it with force'() {
        given:
        client.getModule(_) >> { throw notFound() }
        client.downloadModuleRelease(_, _, _) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        client.getPipelineRelease('nf-core/rnaseq', '2.0.0') >> release('2.0.0', 'new')
        fetcher.fetch(_, _, 'abc', _) >> { args -> writeMain(args[3] as Path, 'x') }
        fetcher.fetch(_, _, 'new', _) >> { args -> writeMain(args[3] as Path, 'new') }
        resolver().installWithDependencies(ref)
        def main = dir().resolve('main.nf')
        main.text = 'edited'
        ModuleChecksum.clearPipeline(dir())

        when:
        resolver().installWithDependencies(ref, '2.0.0')
        then:
        def e = thrown(AbortOperationException)
        e.message.startsWith('Pipeline nf-core/rnaseq has local modifications')
        main.text == 'edited'

        when:
        resolver().installWithDependencies(ref, '2.0.0', true)
        then:
        main.text == 'new'
        ModuleInfo.load(dir(), 'version') == '2.0.0'
        new ModuleStorage(base).getInstalledModule(ref).integrity == ModuleIntegrity.VALID
    }

    def 'replaces an unmodified pipeline with another version'() {
        given:
        client.getModule(_) >> { throw notFound() }
        client.downloadModuleRelease(_, _, _) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        client.getPipelineRelease('nf-core/rnaseq', '2.0.0') >> release('2.0.0', 'new')
        fetcher.fetch(_, _, 'abc', _) >> { args -> writeMain(args[3] as Path, 'old') }
        fetcher.fetch(_, _, 'new', _) >> { args -> writeMain(args[3] as Path, 'new') }
        resolver().installWithDependencies(ref)

        when:
        resolver().installWithDependencies(ref, '2.0.0')
        then:
        dir().resolve('main.nf').text == 'new'
        ModuleInfo.load(dir(), 'version') == '2.0.0'
        Files.list(dir().parent).count() == 1
    }

    def 'a failed install leaves no folder'() {
        given:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))

        when:
        resolver().installWithDependencies(ref)
        then:
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('.git').toFile().mkdirs(); throw new AbortOperationException('boom') }
        thrown(AbortOperationException)
        !Files.exists(dir())
        Files.list(dir().parent).count() == 0

        when:
        resolver().installWithDependencies(ref)
        then:
        1 * fetcher.fetch(*_) >> { args -> writeMain(args[3] as Path, 'x') }
        Files.exists(dir().resolve('main.nf'))
    }

    def 'fails on a revoked release without fetching'() {
        given:
        client.getModule(_) >> { throw notFound() }

        when:
        resolver().installWithDependencies(ref)
        then:
        client.getPipeline(_) >> pipeline(release('1.0.0'), true)
        def e = thrown(AbortOperationException)
        e.message.contains('is revoked')
        0 * fetcher._
        !Files.exists(dir())
    }

    def 'reports a name missing as both a module and a pipeline with the registry detail'() {
        when:
        resolver().installWithDependencies(ref, VERSION)
        then:
        client.getModule(_) >> { throw notFound() }
        client.downloadModuleRelease(_, _, _) >> { throw notFound() }
        client.getPipeline(_) >> { throw new RegistryException('Unable to GET /v1/pipelines/nf-core/rnaseq from any configured registry', new RegistryException('Failed request: x', 404)) }
        def e = thrown(AbortOperationException)
        e.message == MESSAGE

        where:
        VERSION | MESSAGE
        null    | '`nf-core/rnaseq` not found in the registry -- Unable to GET /v1/pipelines/nf-core/rnaseq from any configured registry (Failed request: x)'
        '1.0.0' | '`nf-core/rnaseq@1.0.0` not found in the registry -- Unable to GET /v1/pipelines/nf-core/rnaseq from any configured registry (Failed request: x)'
    }

    def 'fails to install a pipeline without main.nf and leaves no folder'() {
        when:
        resolver().installWithDependencies(ref)
        then:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0')).workflowPath('workflows/rnaseq.nf')
        1 * fetcher.fetch(*_) >> { args -> args[3].resolve('workflows').toFile().mkdirs(); null }
        def e = thrown(AbortOperationException)
        e.message == 'Pipeline nf-core/rnaseq@1.0.0 has no `main.nf` (its main script is `workflows/rnaseq.nf`) and cannot be installed'
        Files.list(dir().parent).count() == 0
    }

    def 'rethrows registry errors other than not found'() {
        when:
        resolver().installWithDependencies(ref)
        then:
        client.getModule(_) >> { throw new RegistryException('Server error', 500) }
        0 * client.getPipeline(_)
        def e = thrown(RegistryException)
        e.statusCode == 500

        when:
        resolver().installWithDependencies(ref)
        then:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> { throw new RegistryException('Unavailable', 503) }
        e = thrown(RegistryException)
        e.statusCode == 503
    }

    def 'a module install with a version makes no extra call'() {
        given:
        def src = Files.createDirectories(base.resolve('src'))
        src.resolve('main.nf').text = 'process RNASEQ { }'
        src.resolve('meta.yml').text = 'name: nf-core/rnaseq\nversion: 1.0.0\n'

        when:
        resolver().installWithDependencies(ref, '1.0.0')
        then:
        1 * client.downloadModuleRelease('nf-core/rnaseq', '1.0.0', _) >> { args -> ModuleStorage.createBundle(src, args[2] as Path); 'oci://x' }
        0 * client._
        0 * fetcher._
        and:
        def installed = new ModuleStorage(base).getInstalledModule(ref)
        !installed.pipeline
        installed.installedVersion == '1.0.0'
    }

    def 'update-deps on a pipeline fails'() {
        given:
        client.getModule(_) >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args -> writeMain(args[3] as Path, 'x') }
        resolver().installWithDependencies(ref)

        when:
        resolver().updateDependencies(ref)
        then:
        def e = thrown(AbortOperationException)
        e.message == 'nf-core/rnaseq is a pipeline -- -update-deps applies to modules only'
    }

    def 'does not vendor the meta.yml dependencies of a pipeline'() {
        when:
        resolver().installWithDependencies(ref)
        then:
        client.getModule('nf-core/rnaseq') >> { throw notFound() }
        client.getPipeline(_) >> pipeline(release('1.0.0'))
        fetcher.fetch(*_) >> { args ->
            args[3].resolve('main.nf').text = 'x'
            args[3].resolve('meta.yml').text = 'name: nf-core/rnaseq\nversion: 1.0.0\nrequires:\n  modules:\n    - nf-core/fastqc@1.0.0\n'
            null
        }
        0 * client.getModule('nf-core/fastqc')
        0 * client.downloadModuleRelease(*_)
        !Files.exists(dir().resolve('modules'))
        new ModuleStorage(base).getInstalledModule(ref).integrity == ModuleIntegrity.VALID
    }
}
