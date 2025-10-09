/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cli

import org.eclipse.jgit.transport.URIish

import java.nio.file.Files

import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import org.eclipse.jgit.api.Git
import spock.lang.IgnoreIf
import spock.lang.Specification

/**
 * Tests for CmdPush command
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class CmdPushTest extends Specification {

    def cleanup() {
        Plugins.stop()
    }

    def 'should fail with no arguments'() {
        given:
        def cmd = new CmdPush()

        when:
        cmd.run()

        then:
        thrown(AbortOperationException)
    }

    def 'should fail when folder does not exist'() {
        given:
        def cmd = new CmdPush(args: ['/nonexistent/folder'], repository: 'https://github.com/test/repo.git')

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Folder does not exist')
    }

    def 'should fail when path is not a directory'() {
        given:
        def tempFile = Files.createTempFile('test', '.txt').toFile()
        def cmd = new CmdPush(args: [tempFile.absolutePath], repository: 'https://github.com/test/repo.git')

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Path is not a directory')

        cleanup:
        tempFile?.delete()
    }

    def 'should fail when no repository specified and no git repo exists'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def cmd = new CmdPush(args: [tempDir.absolutePath])

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('No git repository found')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should normalize repository URLs correctly'() {
        given:
        def cmd = new CmdPush()

        expect:
        cmd.normalizeRepoUrl('https://github.com/user/repo.git') == 'https://github.com/user/repo'
        cmd.normalizeRepoUrl('https://github.com/user/repo') == 'https://github.com/user/repo'
        cmd.normalizeRepoUrl('HTTPS://GITHUB.COM/USER/REPO.GIT') == 'https://github.com/user/repo'
        cmd.normalizeRepoUrl('https://github.com/user/repo/') == 'https://github.com/user/repo'
    }

    def 'should add files to gitignore'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def gitignoreFile = new File(tempDir, '.gitignore')
        def cmd = new CmdPush()

        when:
        cmd.addToGitignore(tempDir, ['file1.txt', 'file2.txt'])

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
        def cmd = new CmdPush()

        when:
        cmd.addToGitignore(tempDir, ['existing.txt', 'new.txt'])

        then:
        def lines = gitignoreFile.readLines()
        lines.count { it == 'existing.txt' } == 1
        lines.contains('new.txt')

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should find work directories'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def workDir = new File(tempDir, 'work')
        workDir.mkdirs()

        def cmd = new CmdPush()

        when:
        def result = cmd.findWorkDirectories(tempDir)

        then:
        result.size() == 1
        result[0] == 'work'

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should not find work directories when none exist'() {
        given:
        def tempDir = Files.createTempDirectory('test').toFile()
        def cmd = new CmdPush()

        when:
        def result = cmd.findWorkDirectories(tempDir)

        then:
        result.isEmpty()

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

        def cmd = new CmdPush(
            args: [tempDir.absolutePath],
            repository: 'https://github.com/correct/repo.git'
        )

        when:
        cmd.run()

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

        def cmd = new CmdPush(
            args: [tempDir.absolutePath],
            repository: 'https://github.com/test/repo.git'
        )

        when:
        cmd.run()

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

        def cmd = new CmdPush(
            args: [tempDir.absolutePath],
            repository: 'https://github.com/test/repo.git',
            revision: 'main'
        )

        when:
        cmd.run()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains("Current branch 'dev' does not match requested branch 'main'")

        cleanup:
        tempDir?.deleteDir()
    }

    def 'should get command name'() {
        given:
        def cmd = new CmdPush()

        expect:
        cmd.getName() == 'push'
    }

    def 'should have default parameter values'() {
        given:
        def cmd = new CmdPush()

        expect:
        cmd.revision == 'main'
        cmd.maxSizeMB == 10
        cmd.message == 'Push from nextflow'
    }
}
