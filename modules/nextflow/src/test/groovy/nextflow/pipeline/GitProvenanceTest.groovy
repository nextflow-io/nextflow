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

import org.eclipse.jgit.api.Git
import spock.lang.Specification
import spock.lang.TempDir

class GitProvenanceTest extends Specification {

    @TempDir
    Path tmp

    private Git repo(String origin = 'https://github.com/org/p.git') {
        final git = Git.init().setDirectory(tmp.toFile()).setInitialBranch('main').call()
        if( origin ) {
            final cfg = git.repository.config
            cfg.setString('remote', 'origin', 'url', origin)
            cfg.save()
        }
        return git
    }

    private void commit(Git git, String file='main.nf', String text='x') {
        final f = tmp.resolve(file)
        Files.createDirectories(f.parent)
        f.text = text
        git.add().addFilepattern(file).call()
        git.commit().setMessage('c').setSign(false).call()
    }

    def 'reads tag at HEAD'() {
        given:
        def git = repo()
        commit(git)
        git.tag().setName('v1.2.0').setAnnotated(false).call()

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.revision == 'v1.2.0'
        p.repositoryUrl == 'https://github.com/org/p'
        p.commitId ==~ /[0-9a-f]{40}/
        p.clean
    }

    def 'uses the branch when no tag points at HEAD'() {
        given:
        def git = repo()
        commit(git)
        git.tag().setName('v1.0.0').setAnnotated(false).call()
        commit(git, 'b.nf')

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.revision == 'main'
    }

    def 'normalizes ssh origin'() {
        given:
        def git = repo('git@github.com:org/p.git')
        commit(git)

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.repositoryUrl == 'https://github.com/org/p'
    }

    def 'not clean with untracked or modified files'() {
        given:
        def git = repo()
        commit(git)

        expect:
        GitProvenance.of(tmp).clean

        when:
        tmp.resolve('new.nf').text = 'y'
        then:
        !GitProvenance.of(tmp).clean

        when:
        tmp.resolve('new.nf').delete()
        tmp.resolve('main.nf').text = 'changed'
        then:
        !GitProvenance.of(tmp).clean

        cleanup:
        git.close()
    }

    def 'clean only considers the pipeline subdirectory'() {
        given:
        def git = repo()
        commit(git, 'sub/main.nf')
        commit(git, 'other/x.nf')
        tmp.resolve('other/x.nf').text = 'changed'

        expect:
        GitProvenance.of(tmp.resolve('sub')).clean
        !GitProvenance.of(tmp).clean

        when:
        tmp.resolve('sub/main.nf').text = 'changed'
        then:
        !GitProvenance.of(tmp.resolve('sub')).clean

        cleanup:
        git.close()
    }

    def 'null outside git'() {
        expect:
        GitProvenance.of(tmp) == null
    }

    def 'null repositoryUrl without origin'() {
        given:
        def git = repo(null)
        commit(git)

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.repositoryUrl == null
        p.revision == 'main'
        p.commitId
        p.clean
    }

    def 'normalizes origin #origin'() {
        given:
        def git = repo(origin)
        commit(git)

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.repositoryUrl == expected
        p.commitId

        where:
        origin                             | expected
        'https://github.com/org/p.git'     | 'https://github.com/org/p'
        'git@github.com:org/p.git'         | 'https://github.com/org/p'
        'ssh://git@host/org/p.git'         | 'https://host/org/p'
        'ssh://git@host:2222/org/p'        | 'https://host/org/p'
        'https://user:tok@host/org/p.git'  | 'https://host/org/p'
        'https://github.com/org/p/'        | 'https://github.com/org/p'
        'https://github.com'               | null
        'file:///tmp/p.git'                | null
        '/tmp/p.git'                       | null
        '../p.git'                         | null
    }

    def 'reads annotated tag at HEAD'() {
        given:
        def git = repo()
        commit(git)
        git.tag().setName('v2.0.0').setAnnotated(true).setMessage('rel').setSigned(false).call()

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.revision == 'v2.0.0'
    }

    def 'uses the commit id on detached HEAD'() {
        given:
        def git = repo()
        commit(git)
        commit(git, 'b.nf')
        git.checkout().setName(git.repository.resolve('HEAD~1').name()).call()

        when:
        def p = GitProvenance.of(tmp)
        git.close()

        then:
        p.revision == p.commitId
    }

    def 'null for a repo without commits'() {
        given:
        def git = repo()

        expect:
        GitProvenance.of(tmp) == null

        cleanup:
        git.close()
    }
}
