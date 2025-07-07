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

package nextflow.scm
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import org.eclipse.jgit.api.Git
import org.eclipse.jgit.lib.Config
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalRepositoryProviderTest extends Specification {

    @Shared
    Path testFolder

    @Shared
    Git repo

    def setup() {
        testFolder = Files.createTempDirectory('test').toAbsolutePath()

        // create the main project
        testFolder.resolve('project_hello').mkdirs()
        def dir = testFolder.resolve('project_hello').toFile()
        def init = Git.init()
        this.repo = init.setDirectory( dir ).call()
        new File(dir, 'main.nf').text = 'main script'
        repo.add().addFilepattern('main.nf').call()
        repo.commit().setSign(false).setMessage('First commit').call()
    }

    def cleanup() {
        testFolder?.deleteDir()
    }

    def 'should return the repo url' () {

        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)
        expect:
        manager.getEndpointUrl() == "file:${testFolder}/project_hello".toString()

    }

    def 'should return content url' () {

        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)
        expect:
        manager.getContentUrl('main.nf') == "file:${testFolder}/project_hello/main.nf".toString()

    }

    def 'should return local clone url'() {
        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)
        expect:
        manager.getCloneUrl() == "file:${testFolder}/project_hello/.git".toString()
    }

    def 'should return local clone url for bare repo'() {

        given:
        def cloneDir = new File("$testFolder/bare_repo")
        Git
                .cloneRepository()
                .setBare(true)
                .setURI("file:${testFolder}/project_hello")
                .setDirectory(cloneDir)
                .call()

        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('bare_repo', config)

        expect:
        manager.getCloneUrl() == "file:${testFolder}/bare_repo".toString()
    }

    def 'should read file bytes' () {
        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)
        expect:
        new String(manager.readBytes('main.nf')) == "main script"
    }

    def 'should read file from bare repository' () {

        given:
        println testFolder
        def cloneDir = new File("$testFolder/bare_repo")
        Git
                .cloneRepository()
                .setBare(true)
                .setURI("file:${testFolder}/project_hello")
                .setDirectory(cloneDir)
                .call()

        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('bare_repo', config)
        expect:
        new String(manager.readBytes('main.nf')) == "main script"

    }

    def 'should list tags' () {
        given:
        def dir = testFolder.toFile()
        def ref1 = repo.tag().setName('tag_1').setMessage('First tag').call()
        // add new file
        new File(dir, 'foo.nf').text = 'foo script'
        repo.add().addFilepattern('foo.nf').call()
        repo.commit().setSign(false).setMessage('Second commit').call()
        def ref2 = repo.tag().setName('tag_2').setMessage('Second tag').call()

        and:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)

        when:
        def tags = manager.getTags()
        then:
        tags.size() == 2
        and:
        tags[0].name == 'tag_1'
        tags[0].getCommitId() == ref1.getObjectId().name()
        and:
        tags[1].name == 'tag_2'
        tags[1].getCommitId() == ref2.getObjectId().name()
    }

    def 'should list branches' () {
        given:
        def dir = testFolder.toFile()
        def ref1 = repo.branchCreate().setName('branch_1').call()
        // add new file
        new File(dir, 'foo.nf').text = 'foo script'
        repo.add().addFilepattern('foo.nf').call()
        repo.commit().setSign(false).setMessage('Second commit').call()
        def ref2 = repo.branchCreate().setName('branch_2').call()

        // Use this user's custom defaultBranch name if set in ~/.gitconfig
        def defaultBranch = 'master'
        def gitconfig = Paths.get(System.getProperty('user.home'),'.gitconfig');
        if(gitconfig.exists()) {
            def config = new Config()
            config.fromText(gitconfig.text)
            defaultBranch = config.getString('init', null, 'defaultBranch') ?: 'master'
        }

        and:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('project_hello', config)

        when:
        def branches = manager.getBranches()
        then:
        branches.size() == 3
        and:
        branches.find { it.name == defaultBranch }
        and:
        branches.find { it.name == 'branch_1' }.commitId == ref1.getObjectId().name()
        and:
        branches.find { it.name == 'branch_2' }.commitId == ref2.getObjectId().name()
    }
}
