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
        def param = [:] as FileOutParam

        expect:
        param.relativize( Paths.get('x'), workDir ) == 'x'
        param.relativize( Paths.get('hello.txt'), workDir ) == 'hello.txt'
        param.relativize( Paths.get('sub/dir/hello.txt'), workDir ) == 'sub/dir/hello.txt'
        param.relativize( Paths.get('/a/b/c/x'), workDir ) == 'x'
        param.relativize( Paths.get('/a/b/c/hello.txt'), workDir ) == 'hello.txt'
        param.relativize( Paths.get('/a/b/c/some/path/hello.txt'), workDir ) == 'some/path/hello.txt'

        when:
        param.relativize( Paths.get('/c/b/a/hello.txt'), workDir )
        then:
        thrown(IllegalFileException)

        when:
        param.relativize( Paths.get('/a/b/c'), workDir )
        then:
        thrown(IllegalFileException)
    }

    def 'should escape glob' () {

        given:
        def workDir = Paths.get('/a/b/c')
        def param1 = [:] as FileOutParam
        def param2 = [:] as FileOutParam; param2.glob = false
        def param3 = [:] as FileOutParam; param3.glob = true

        expect:
        param1.relativize( Paths.get('hello[a-b].txt'), workDir ) == 'hello\\[a-b\\].txt'
        param1.relativize( Paths.get('/a/b/c/some/path/hello[a-b].txt'), workDir ) == 'some/path/hello\\[a-b\\].txt'

        param2.relativize( Paths.get('hello[a-b].txt'), workDir ) == 'hello[a-b].txt'
        param2.relativize( Paths.get('/a/b/c/some/path/hello[a-b].txt'), workDir ) == 'some/path/hello[a-b].txt'

        param3.relativize( Paths.get('hello[a-b].txt'), workDir ) == 'hello\\[a-b\\].txt'
        param3.relativize( Paths.get('/a/b/c/some/path/hello[a-b].txt'), workDir ) == 'some/path/hello\\[a-b\\].txt'

    }


    def 'should return a name relative to the workDir (with string)' () {

        given:
        def param = [:] as FileOutParam
        def workDir = Paths.get('/a/b/c/')

        expect:
        param.relativize( 'x', workDir ) == 'x'
        param.relativize( 'hello.txt', workDir ) == 'hello.txt'
        param.relativize( 'sub/dir/hello.txt', workDir ) == 'sub/dir/hello.txt'
        param.relativize( '/a/b/c/x', workDir ) == 'x'
        param.relativize( '/a/b/c/hello.txt', workDir ) == 'hello.txt'
        param.relativize( '/a/b/c/some/path/hello.txt', workDir ) == 'some/path/hello.txt'

        when:
        param.relativize( '/c/b/a/hello.txt', workDir )
        then:
        thrown(IllegalFileException)

        when:
        param.relativize( '/a/b/c', workDir )
        then:
        thrown(IllegalFileException)

        when:
        param.relativize( '/a/b/c/', workDir )
        then:
        thrown(IllegalFileException)

    }


}
