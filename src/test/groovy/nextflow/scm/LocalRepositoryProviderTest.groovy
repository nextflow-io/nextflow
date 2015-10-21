/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

    def setupSpec() {
        testFolder = Files.createTempDirectory('test').toAbsolutePath()

        // create the main project
        testFolder.resolve('test/pipe').mkdirs()
        def dir = testFolder.resolve('test/pipe').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        new File(dir, 'main.nf').text = 'main script'
        repo.add().addFilepattern('*').call()
        repo.commit().setAll(true).setMessage('First commit').call()

    }

    def cleanupSpec() {
        testFolder?.deleteDir()
    }

    def 'should return the repo url' () {

        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('nextflow-io/hello', config)
        expect:
        manager.getEndpointUrl() == "file:${testFolder}/nextflow-io/hello".toString()

    }

    def 'should return content url' () {

        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('nextflow-io/hello', config)
        expect:
        manager.getContentUrl('main.nf') == "file:${testFolder}/nextflow-io/hello/main.nf".toString()

    }

    def 'should return local clone url'() {
        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('nextflow-io/hello', config)
        expect:
        manager.getCloneUrl() == "file:${testFolder}/nextflow-io/hello/.git".toString()
    }

    def 'should read file bytes' () {
        given:
        def config = new ProviderConfig('local', [path: testFolder])
        def manager = new LocalRepositoryProvider('test/pipe', config)
        expect:
        new String(manager.readBytes('main.nf')) == "main script"
    }


}
