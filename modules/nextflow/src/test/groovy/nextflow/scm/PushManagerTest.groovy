/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.scm

import nextflow.exception.AbortOperationException
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.transport.URIish
import spock.lang.Specification

import java.nio.file.Files

/**
 * Tests for PushManager
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class PushManagerTest extends Specification {

    def 'should normalize repository URLs correctly'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def manager = new PushManager(tempDir, false, 10, false, false)

        expect:
        manager.normalizeRepoUrl('https://github.com/user/repo.git') == 'https://github.com/user/repo'
        manager.normalizeRepoUrl('https://github.com/user/repo') == 'https://github.com/user/repo'
        manager.normalizeRepoUrl('HTTPS://GITHUB.COM/USER/REPO.GIT') == 'https://github.com/user/repo'
        manager.normalizeRepoUrl('https://github.com/user/repo/') == 'https://github.com/user/repo'

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should add files to gitignore'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def gitignoreFile = new File(tempDir, '.gitignore')
        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.addToGitignore(['file1.txt', 'file2.txt'])

        then:
        gitignoreFile.exists()
        def content = gitignoreFile.text
        content.contains('file1.txt')
        content.contains('file2.txt')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should not duplicate entries in gitignore'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def gitignoreFile = new File(tempDir, '.gitignore')
        gitignoreFile.text = 'existing.txt\n'
        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.addToGitignore(['existing.txt', 'new.txt'])

        then:
        def lines = gitignoreFile.readLines()
        lines.count { it == 'existing.txt' } == 1
        lines.contains('new.txt')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should detect existing git repository'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        when:
        def manager = new PushManager(tempDir, false, 10, false, false)
        def result1 = manager.isLocalGit()

        then:
        !result1

        when:
        Git.init().setDirectory(tempDir).call().close()
        def result2 = manager.isLocalGit()

        then:
        result2

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when existing repo has wrong remote'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize a git repo with a different remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/wrong/repo.git'))
            .call()
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.push('https://github.com/correct/repo.git', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Repository URL not found in remotes')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when repo is in detached HEAD state'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize repo and create a commit
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create a test file and commit
        new File(tempDir, 'test.txt').text = 'test content'
        git.add().addFilepattern('.').call()
        def commit = git.commit().setMessage('initial commit').call()

        // Checkout to detached HEAD
        git.checkout().setName(commit.name()).call()
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('detached HEAD state')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when current branch does not match requested branch'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize repo
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create initial commit on main branch
        new File(tempDir, 'test.txt').text = 'test content'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()

        // Create and checkout dev branch
        git.checkout().setCreateBranch(true).setName('dev').call()
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.push('https://github.com/test/repo.git', 'main')

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Current branch 'dev' does not match requested branch 'main'")

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when no git repository found and commit is false'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No git repository found')
        e.message.contains('commit')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should resolve repository from single remote'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with single remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        def repo = manager.resolveRepository()

        then:
        repo == 'https://github.com/test/repo.git'

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail to resolve repository when no git repo exists'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.resolveRepository()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No git repository found')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail to resolve repository when no remotes configured'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo without remotes
        Git.init().setDirectory(tempDir).call().close()

        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.resolveRepository()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No remotes configured')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when nextflow.config is untracked and commit is false'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create an initial commit to avoid empty repo issues
        new File(tempDir, 'README.md').text = 'test'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()

        // Create untracked nextflow.config file
        new File(tempDir, 'nextflow.config').text = 'process.executor = "local"'
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, true)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('nextflow.config')
        e.message.contains('untracked')
        e.message.contains('commit')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should not fail when nextflow.config is untracked but commit is true'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create an initial commit
        new File(tempDir, 'README.md').text = 'test'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()

        // Create untracked nextflow.config file
        new File(tempDir, 'nextflow.config').text = 'process.executor = "local"'
        git.close()

        def manager = new PushManager(tempDir, true, 10, false,false)

        // Note: This test verifies that checkNextflowConfig doesn't throw an exception
        // The push operation itself will fail because AssetManager.upload is not mocked,
        // but we're only testing that checkNextflowConfig passes validation
        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        // The exception thrown should not be about nextflow.config being untracked
        def e = thrown(Exception)
        !e.message?.contains('nextflow.config')
        !e.message?.contains('untracked')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should fail when nextflow.config is ignored and not in index'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create initial commit
        new File(tempDir, 'README.md').text = 'test'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()

        // Create .gitignore and ignore nextflow.config
        new File(tempDir, '.gitignore').text = 'nextflow.config\n'
        git.add().addFilepattern('.gitignore').call()
        git.commit().setMessage('add gitignore').call()

        // Create nextflow.config file (it will be ignored)
        new File(tempDir, 'nextflow.config').text = 'process.executor = "local"'
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, true)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('nextflow.config')
        e.message.contains('ignored')
        e.message.contains('.gitignore')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should not fail when nextflow.config is properly tracked'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create and commit nextflow.config
        new File(tempDir, 'nextflow.config').text = 'process.executor = "local"'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()
        git.close()

        def manager = new PushManager(tempDir, false, 10, false, true)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        // The exception thrown should not be about nextflow.config
        // It will fail at the AssetManager.upload stage since we're not mocking
        def e = thrown(Exception)
        !e.message?.contains('nextflow.config')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should not fail when checkConfigExist is false and nextflow.config is untracked'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()

        // Initialize git repo with remote
        def git = Git.init().setDirectory(tempDir).call()
        git.remoteAdd()
            .setName('origin')
            .setUri(new URIish('https://github.com/test/repo.git'))
            .call()

        // Create an initial commit
        new File(tempDir, 'README.md').text = 'test'
        git.add().addFilepattern('.').call()
        git.commit().setMessage('initial commit').call()

        // Create untracked nextflow.config file
        new File(tempDir, 'nextflow.config').text = 'process.executor = "local"'
        git.close()

        // Create manager with checkConfigExist = false
        def manager = new PushManager(tempDir, false, 10, false, false)

        when:
        manager.push('https://github.com/test/repo.git', null)

        then:
        // The exception thrown should not be about nextflow.config validation
        def e = thrown(Exception)
        !e.message?.contains('nextflow.config')
        !e.message?.contains('untracked')

        cleanup:
        tempDir?.deleteDir()
    }
}
