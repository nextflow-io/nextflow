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

package nextflow
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.util.ArrayTuple
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowTest extends Specification {


    def testFile() {

        expect:
        Nextflow.file('file.log').toFile() == new File('file.log').canonicalFile
        Nextflow.file('relative/file.test').toFile() == new File( new File('.').canonicalFile, 'relative/file.test')
        Nextflow.file('/user/home/file.log').toFile() == new File('/user/home/file.log')
        Nextflow.file('~').toFile() == new File( System.getProperty('user.home') )
        Nextflow.file('~/file.test').toFile() == new File( System.getProperty('user.home'), 'file.test' )
        Nextflow.file('~file.name').toFile() == new File('~file.name').canonicalFile

    }



    def testFile2() {

        when:
        def current = new File('.').canonicalPath

        then:
        Nextflow.file('hola').toString() == current + '/hola'
        Nextflow.file( new File('path/file.txt') ).toString() == current + '/path/file.txt'
        Nextflow.file( Paths.get('some/path') ).toString() == current + '/some/path'
        Nextflow.file( '/abs/path/file.txt' ) == Paths.get('/abs/path/file.txt')


    }

    def testFileWithWildcards() {

        setup:
        def folder = Files.createTempDirectory('nxftest') .toAbsolutePath()
        folder.resolve('hola1').text = 'abc'
        folder.resolve('helo2').text = 'abc'
        folder.resolve('hello3').text = 'abc'

        expect:
        Nextflow.file(folder.toString() + '/ciao*') == []
        Nextflow.file(folder.toString() + '/hel*').sort().collect { it.name } == [ 'helo2','hello3' ].sort()
        Nextflow.file(folder.toString() + '/hol??').collect { it.name } == [ 'hola1' ]


        cleanup:
        folder?.delete()

    }


    def testTuple() {

        expect:
        Nextflow.tuple(1) == new ArrayTuple([1])
        Nextflow.tuple(1,2,3) == new ArrayTuple([1,2,3])
        Nextflow.tuple([4,5,6]) == new ArrayTuple([4,5,6])
        Nextflow.tuple([4,5,6], [1,2]) == new ArrayTuple([[4,5,6], [1,2]])
    }


}
