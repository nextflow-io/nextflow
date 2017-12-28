/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
