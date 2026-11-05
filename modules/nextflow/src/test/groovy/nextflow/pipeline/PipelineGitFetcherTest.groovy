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

import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import spock.lang.Specification
import spock.lang.TempDir

class PipelineGitFetcherTest extends Specification {

    @TempDir
    Path dir

    String sha1
    String sha2
    String url

    def setup() {
        def repo = dir.resolve('repo')
        def git = Git.init().setDirectory(repo.toFile()).call()
        repo.resolve('main.nf').text = 'one'
        git.add().addFilepattern('.').call()
        sha1 = git.commit().setMessage('first').setSign(false).call().name()
        git.tag().setName('1.0').call()
        repo.resolve('main.nf').text = 'two'
        git.add().addFilepattern('.').call()
        sha2 = git.commit().setMessage('second').setSign(false).call().name()
        git.close()
        url = repo.toUri().toString()
    }

    def 'fetches the tag and removes .git'() {
        given:
        def target = dir.resolve('target')
        when:
        new PipelineGitFetcher().fetch(url, '1.0', sha1, target)
        then:
        target.resolve('main.nf').text == 'one'
        !Files.exists(target.resolve('.git'))
    }

    def 'fetches a commit sha revision'() {
        given:
        def target = dir.resolve('target')
        when:
        new PipelineGitFetcher().fetch(url, sha2, sha2, target)
        then:
        target.resolve('main.nf').text == 'two'
        !Files.exists(target.resolve('.git'))
    }

    def 'fetches the pinned commit when the revision is blank'() {
        given:
        def target = dir.resolve('target')
        when:
        new PipelineGitFetcher().fetch(url, revision, sha1, target)
        then:
        target.resolve('main.nf').text == 'one'
        !Files.exists(target.resolve('.git'))
        where:
        revision << [null, '']
    }

    def 'fails when the revision no longer points to the pinned commit'() {
        when:
        new PipelineGitFetcher().fetch(url, '1.0', sha2, dir.resolve('target'))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('does not match the pinned commit')
    }

    def 'fails on an unknown revision'() {
        when:
        new PipelineGitFetcher().fetch(url, 'no-such-tag', sha1, dir.resolve('target'))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Cannot fetch revision')
    }
}
