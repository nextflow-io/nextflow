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

import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StagePathTest extends Specification {

    def testToString() {

        when:
        def real = Files.createTempDirectory('test')
        def path = new StagePath(real)

        then:
        path.toString() == real.getFileName().toString()

        cleanup:
        real.deleteDir()

    }


    def testTextProperty() {

        given:
        def folder = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'Hello'
        then:
        stage.text == 'Hello'

        when:
        stage.text = 'Ciao!'
        then:
        real.text == 'Ciao!'

        cleanup:
        folder.deleteDir()

    }

    def testReadLines() {

        given:
        def folder = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'a\nbb\nccc'
        then:
        stage.readLines() == ['a','bb','ccc']

        cleanup:
        folder.deleteDir()

    }

    def testCopyTo() {

        given:
        def folder = Files.createTempDirectory('test')
        def target = Files.createTempDirectory('test')
        def real = folder.resolve('name')
        def stage = new StagePath(real)

        when:
        real.text = 'a\nbb\nccc'
        stage.copyTo(target)
        then:
        target.resolve('name').text == 'a\nbb\nccc'

        cleanup:
        folder.deleteDir()
        target.deleteDir()

    }
}
