/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow.script

import java.nio.file.Paths

import nextflow.exception.IllegalFileException
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileOutParamTest extends Specification {

    def 'should return a name relative to the workDir' () {

        given:
        def workDir = Paths.get('/a/b/c')

        expect:
        FileOutParam.relativize( Paths.get('x'), workDir ) == 'x'
        FileOutParam.relativize( Paths.get('hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativize( Paths.get('sub/dir/hello.txt'), workDir ) == 'sub/dir/hello.txt'
        FileOutParam.relativize( Paths.get('/a/b/c/x'), workDir ) == 'x'
        FileOutParam.relativize( Paths.get('/a/b/c/hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativize( Paths.get('/a/b/c/some/path/hello.txt'), workDir ) == 'some/path/hello.txt'

        when:
        FileOutParam.relativize( Paths.get('/c/b/a/hello.txt'), workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( Paths.get('/a/b/c'), workDir )
        then:
        thrown(IllegalFileException)
    }


    def 'should return a name relative to the workDir (with string)' () {

        given:
        def workDir = Paths.get('/a/b/c/')

        expect:
        FileOutParam.relativize( 'x', workDir ) == 'x'
        FileOutParam.relativize( 'hello.txt', workDir ) == 'hello.txt'
        FileOutParam.relativize( 'sub/dir/hello.txt', workDir ) == 'sub/dir/hello.txt'
        FileOutParam.relativize( '/a/b/c/x', workDir ) == 'x'
        FileOutParam.relativize( '/a/b/c/hello.txt', workDir ) == 'hello.txt'
        FileOutParam.relativize( '/a/b/c/some/path/hello.txt', workDir ) == 'some/path/hello.txt'

        when:
        FileOutParam.relativize( '/c/b/a/hello.txt', workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( '/a/b/c', workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( '/a/b/c/', workDir )
        then:
        thrown(IllegalFileException)

    }


}
