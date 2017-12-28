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

package nextflow.file

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Session
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilePorterTest extends Specification {

    def 'should get the max threads value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxThreads() == Runtime.getRuntime().availableProcessors()

        when:
        new Session([filePorter: [maxThreads: 99]])
        then:
        FilePorter.getMaxThreads() == 99

    }

    def 'should get the max retries value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxRetries() == 3

        when:
        new Session([filePorter: [maxRetries: 88]])
        then:
        FilePorter.getMaxRetries() == 88

    }

    def 'should copy foreign files' () {

        given:
        new Session()
        def folder = Files.createTempDirectory('test')
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')
        def local = Paths.get('local.txt')
        def files = [foo: local, bar: foreign1, baz: foreign2]

        when:
        def d = new FilePorter(folder)
        def result = d.stageForeignFiles(files)
        then:
        result.foo ==  Paths.get('local.txt')

        result.bar.name == 'hola.txt'
        result.bar.text == 'hola mundo!'
        result.bar.fileSystem == FileSystems.default

        result.baz.name == 'ciao.txt'
        result.baz.text == 'ciao mondo!'
        result.baz.fileSystem == FileSystems.default

        cleanup:
        folder?.deleteDir()
    }
}
