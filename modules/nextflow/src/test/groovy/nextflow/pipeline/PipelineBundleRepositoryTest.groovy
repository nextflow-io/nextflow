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

import io.seqera.npr.client.RegistryClient
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Constants
import spock.lang.Specification
import spock.lang.TempDir

class PipelineBundleRepositoryTest extends Specification {

    @TempDir
    Path folder

    /** A bundle with a main.nf printing the version and an executable script */
    private Map bundle(String version) {
        final src = Files.createDirectories(folder.resolve("src-$version"))
        src.resolve('main.nf').text = "println '$version'"
        final bin = Files.createDirectories(src.resolve('bin')).resolve('run.sh')
        bin.text = 'echo'
        bin.toFile().setExecutable(true)
        final target = folder.resolve("b-${version}.tar.gz")
        final checksum = PipelineBundle.create(src, target)
        return [path: target, checksum: checksum]
    }

    def 'should add each bundle release as a tagged commit, with the default branch on the latest'() {
        given:
        def repo = folder.resolve('bundles.git')
        def client = Mock(RegistryClient)
        def bundles = ['1.0.0', '2.0.0', '1.5.0'].collectEntries { [(it): bundle(it)] }
        def release = { String v -> new PipelineInfo(name: 'acme/hello', version: v, artifactChecksum: bundles[v].checksum as String) }

        when:
        ['1.0.0', '2.0.0', '1.5.0'].each { PipelineBundleRepository.add(repo, client, release(it)) }

        then:
        3 * client.downloadPipelineRelease('acme/hello', _, _) >> { args -> Files.copy(bundles[args[1]].path as Path, args[2] as Path); null }

        and:
        def git = Git.open(repo.toFile())
        git.tagList().call()*.name == ['refs/tags/1.0.0', 'refs/tags/1.5.0', 'refs/tags/2.0.0']
        git.repository.resolve('refs/heads/master') == git.repository.resolve('refs/tags/2.0.0^{commit}')
        new String(git.repository.open(git.repository.resolve('1.5.0:main.nf')).bytes) == "println '1.5.0'"
        git.repository.resolve('1.0.0:bin/run.sh') != null
        // every version is reachable from the releases branch, so a fetch of any branch gets all tags
        def releases = git.log().add(git.repository.resolve('refs/heads/releases')).call()*.name
        ['1.0.0', '1.5.0', '2.0.0'].every { git.repository.resolve("refs/tags/${it}^{commit}").name in releases }

        when: 'adding a version already in the repository'
        PipelineBundleRepository.add(repo, client, release('1.0.0'))

        then:
        0 * client.downloadPipelineRelease(*_)

        cleanup:
        git?.close()
    }

    def 'should order versions as semver'() {
        expect:
        Math.signum(PipelineBundleRepository.compare(A, B)) == EXPECTED

        where:
        A            | B            | EXPECTED
        '2.0.0'      | '1.10.0'     | 1
        '1.2.0'      | '1.10.0'     | -1
        '2.0.0-rc1'  | '2.0.0'      | -1
        '2.0.0-rc2'  | '2.0.0-rc1'  | 1
        '1.0.0'      | '1.0.0'      | 0
    }

    def 'should refuse a bundle that does not match the registry checksum'() {
        given:
        def repo = folder.resolve('bundles.git')
        def client = Mock(RegistryClient)
        def b = bundle('1.0.0')

        when:
        PipelineBundleRepository.add(repo, client, new PipelineInfo(name: 'acme/hello', version: '1.0.0', artifactChecksum: 'sha256:00'))

        then:
        1 * client.downloadPipelineRelease(*_) >> { args -> Files.copy(b.path as Path, args[2] as Path); null }
        thrown(nextflow.exception.AbortOperationException)
        Git.open(repo.toFile()).withCloseable { it.repository.findRef(Constants.R_TAGS + '1.0.0') } == null
    }
}
