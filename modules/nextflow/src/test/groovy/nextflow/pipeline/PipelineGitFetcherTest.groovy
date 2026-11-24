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

    // the commits are created in setup, so the tables name them
    private String sha(String value) {
        value == 'SHA1' ? sha1 : value == 'SHA2' ? sha2 : value
    }

    def 'fetches the pinned commit, else the revision, and removes .git'() {
        given:
        def target = dir.resolve('target')
        when:
        def commit = new PipelineGitFetcher().fetch(url, sha(REVISION), sha(COMMIT_ID), target)
        then:
        commit == sha(EXPECTED)
        target.resolve('main.nf').text == (EXPECTED == 'SHA1' ? 'one' : 'two')
        !Files.exists(target.resolve('.git'))
        where:
        REVISION | COMMIT_ID | EXPECTED
        '1.0'    | 'SHA1'    | 'SHA1'   // a tag
        'SHA2'   | 'SHA2'    | 'SHA2'   // a commit sha
        null     | 'SHA1'    | 'SHA1'   // no revision
        ''       | 'SHA1'    | 'SHA1'
        'master' | 'SHA1'    | 'SHA1'   // the revision has moved
        '1.0'    | null      | 'SHA1'   // no commit id
        '1.0'    | ''        | 'SHA1'
    }

    def 'fails clearly when there is neither a revision nor a commit id'() {
        when:
        new PipelineGitFetcher().fetch(url, null, commitId, dir.resolve('target'))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('neither a git revision nor a commit id')
        where:
        commitId << [null, '']
    }

    def 'fails on an unknown revision'() {
        when:
        new PipelineGitFetcher().fetch(url, 'no-such-tag', null, dir.resolve('target'))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Cannot fetch revision 'no-such-tag'")
    }

    def 'fails on an unknown commit'() {
        when:
        new PipelineGitFetcher().fetch(url, '1.0', 'f' * 40, dir.resolve('target'))
        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Cannot fetch commit")
    }
}
