/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import org.eclipse.jgit.api.Git
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalRepositoryProviderTest extends Specification {

    @Shared
    Path testFolder

    def setup() {
        testFolder = Files.createTempDirectory('test').toAbsolutePath()

        // create the main project
        testFolder.resolve('project_hello').mkdirs()
        def dir = testFolder.resolve('project_hello').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        new File(dir, 'main.nf').text = 'main script'
        repo.add().addFilepattern('main.nf').call()
        repo.commit().setMessage('First commit').call()

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


}
