/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
