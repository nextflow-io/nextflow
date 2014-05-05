/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.extension

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Channel
import nextflow.Session
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowExtensions2Test extends Specification {

    @Shared
    def Session session

    def setupSpec() {
        session = new Session()
        session.workDir = Files.createTempDirectory('hello')
    }

    def cleanupSpec() {
        session.workDir.deleteDir()
    }

    def testAppendFileString() {

        when:
        def result = Channel
                .from('alpha','beta','gamma')
                .appendFile { it == 'beta' ? ['file2', it.reverse() ] : ['file1',it] }
                .toSortedList { it.name }

        List<Path> list = result.val

        then:
        list[0].name == 'file1'
        list[0].text == 'alphagamma'

        list[1].name == 'file2'
        list[1].text == 'ateb'

    }


    def testAppendFileWithFiles() {


        given:
        def file1 = Files.createTempDirectory('temp').resolve('A')
        file1.text = 'alpha\nbeta'

        def file2 = Files.createTempDirectory('temp').resolve('B')
        file2.text = 'Hello\nworld'

        def file3 = Files.createTempDirectory('temp').resolve('A')
        file3.text = 'xyz'

        when:
        def list = Channel
                .from(file1,file2,file3)
                .appendFile()
                .toSortedList { it.name }
                .getVal() as List<Path>

        then:
        list[0].name == 'A'
        list[0].text == 'alpha\nbetaxyz'

        list[1].name == 'B'
        list[1].text == 'Hello\nworld'


        when:
        list = Channel
                .from(file1,file2,file3)
                .appendFile(newLine:true)
                .toSortedList { it.name }
                .getVal() as List<Path>


        then:
        list[0].name == 'A'
        list[0].text == 'alpha\nbeta\nxyz\n'

        list[1].name == 'B'
        list[1].text == 'Hello\nworld\n'

    }

}
