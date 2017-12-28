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

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilePatternSplitterTest extends Specification {

    @Unroll
    def 'test isGlobPattern #str' () {

        when:
        def parser = FilePatternSplitter.glob().parse(str)
        then:
        parser.isPattern() == result
        parser.fileName == pattern

        where:
        str             | result        | pattern
        'hola'          | false         | 'hola'
        '1-2-3'         | false         | '1-2-3'
        'hello.txt'     | false         | 'hello.txt'
        'hello{x'       | false         | 'hello{x'
        'hello\\{x'     | false         | 'hello\\{x'
        'hello[x'       | false         | 'hello[x'
        'hello(a)'      | false         | 'hello(a)'
        'some/path'     | false         | 'path'
        '*'             | true          | '*'
        'hola*'         | true          | 'hola*'
        'hola?'         | true          | 'hola?'
        '?'             | true          | '?'
        'hola[a]'       | true          | 'hola[a]'
        'hola[a-z]'     | true          | 'hola[a-z]'
        'hola{a,b}'     | true          | 'hola{a,b}'
        'hola{a,}'      | true          | 'hola{a,}'
        'hola{,a}'      | true          | 'hola{,a}'
        'hola{,}'       | true          | 'hola{,}'
        'hola{a}'       | false         | 'hola{a}'
        'hola[]'        | false         | 'hola[]'
        '\\*'           | false         | '\\*'
        'hola\\*'       | false         | 'hola\\*'
        'hola\\\\*'     | false         | 'hola\\\\*'
        'hola\\?'       | false         | 'hola\\?'
        '\\?'           | false         | '\\?'
        'hola\\[a\\]'   | false         | 'hola\\[a\\]'
        'hola\\[a-z\\]' | false         | 'hola\\[a-z\\]'
        'hola\\{a,b\\}' | false         | 'hola\\{a,b\\}'
        'hola{a,b\\}'   | false         | 'hola{a,b\\}'
        'hola\\{a,b}'   | false         | 'hola\\{a,b}'
        'hola\\a'       | false         | 'hola\\a'
        'hola\\\\a'     | false         | 'hola\\\\a'
        '/tmp/*_{1,2}.fa'   | true      | '*_{1,2}.fa'
        '/tmp/\\*_{1,2}.fa' | true      | '\\*_{1,2}.fa'
    }

    @Unroll
    def 'should split path components for: #str' () {
        given:
        def parser = FilePatternSplitter.glob()

        when:
        parser.parse(str)
        then:
        parser.parent == folder
        parser.fileName == pattern
        parser.scheme == scheme

        where:
        str                                     | folder            | pattern               | scheme
        '/some/file/name.txt'                   | '/some/file/'     | 'name.txt'            | null
        '/some/file/na*.txt'                    | '/some/file/'     | 'na*.txt'             | null
        '/some/file/na??.txt'                   | '/some/file/'     | 'na??.txt'            | null
        '/some/file/*.txt'                      | '/some/file/'     | '*.txt'               | null
        '/some/file/?.txt'                      | '/some/file/'     | '?.txt'               | null
        '/some/file/*'                          | '/some/file/'     | '*'                   | null
        '/some/file/'                           | '/some/file/'     | ''                    | null
        'path/filename.txt'                     | 'path/'           | 'filename.txt'        | null
        'filename.txt'                          | './'              | 'filename.txt'        | null
        './file.txt'                            | './'              | 'file.txt'            | null
        '/some/file/**/*.txt'                   | '/some/file/'     | '**/*.txt'            | null
        'dxfs:///some/file/**/*.txt'            | '/some/file/'     | '**/*.txt'            | 'dxfs'
        'dxfs://some/file/**/*.txt'             | 'some/file/'      | '**/*.txt'            | 'dxfs'
        'dxfs://*.txt'                          | './'              | '*.txt'               | 'dxfs'
        'dxfs:///*.txt'                         | '/'               | '*.txt'               | 'dxfs'
        'dxfs:///**/*.txt'                      | '/'               | '**/*.txt'            | 'dxfs'
        'file{a,b}'                             | './'              | 'file{a,b}'           | null
        'test/data/file{a,b}'                   | 'test/data/'      | 'file{a,b}'           | null
        'test/{file1,file2}'                    | 'test/'           | '{file1,file2}'       | null
        '{file1,file2}'                         | './'              | '{file1,file2}'       | null
        '{test/file1,data/file2}'               | './'              | '{test/file1,data/file2}' | null
        'data/{p/file1,q/file2}'                | 'data/'           | '{p/file1,q/file2}'   | null
        'data/{p,q}/file'                       | 'data/'           | '{p,q}/file'          | null
        'test/data/file[a-b]'                   | 'test/data/'      | 'file[a-b]'           | null
        'test/data[a-b]/file'                   | 'test/'           | 'data[a-b]/file'      | null
        '/some/path\\[a-b\\]/data{a,b}/file\\?' | '/some/path[a-b]/'| 'data{a,b}/file\\?'   | null

    }

    def 'should strip glob escape chars' () {

        expect:
        FilePatternSplitter.glob().strip('hello')  == 'hello'
        FilePatternSplitter.glob().strip('hello-world')  == 'hello-world'
        FilePatternSplitter.glob().strip('hello\\world')  == 'hello\\world'
        FilePatternSplitter.glob().strip('hello/world')  == 'hello/world'
        FilePatternSplitter.glob().strip('hello world?')  == 'hello world?'
        FilePatternSplitter.glob().strip('[hello world]')  == '[hello world]'
        FilePatternSplitter.glob().strip('\\[hello world\\]')  == '[hello world]'
        FilePatternSplitter.glob().strip('hello\\[-\\]world')  == 'hello[-]world'
        FilePatternSplitter.glob().strip('{hello,world}')  == '{hello,world}'
        FilePatternSplitter.glob().strip('\\{hello,world\\}')  == '{hello,world}'
        FilePatternSplitter.glob().strip('hello\\{,\\}world')  == 'hello{,}world'
        FilePatternSplitter.glob().strip('file-\\?.txt')  == 'file-?.txt'
        FilePatternSplitter.glob().strip('file-\\*.txt')  == 'file-*.txt'
        FilePatternSplitter.glob().strip('file-\\\\*.txt')  == 'file-\\*.txt'
    }


    def 'should escape glob characters' () {
        expect:
        FilePatternSplitter.glob().escape('hello') == 'hello'
        FilePatternSplitter.glob().escape('/some/path') == '/some/path'
        FilePatternSplitter.glob().escape('/some/file-*') == '/some/file-\\*'
        FilePatternSplitter.glob().escape('/some/*.*') == '/some/\\*.\\*'
        FilePatternSplitter.glob().escape('/some/file.?') == '/some/file.\\?'
        FilePatternSplitter.glob().escape('/some/file[a-b]') == '/some/file\\[a-b\\]'
        FilePatternSplitter.glob().escape('/some/file{a,b}') == '/some/file\\{a,b\\}'
    }


    @Unroll
    def 'should find first meta char for #str' () {

        given:
        def parser = FilePatternSplitter.glob()

        expect:
        parser.firstMetaIndex(str) == index

        where:
        str             | index
        'hello'         | -1
        'hell*'         | 4
        'hell?-*'       | 4
        '*ello??'       | 0
    }

    @Unroll
    def 'should convert the path to canonical form #str' () {

        given:
        def parser = FilePatternSplitter.glob()

        expect:
        parser.replaceMetaChars(str,'X' as char) == expected

        where:
        str             | expected
        'hello'         | 'hello'
        'hell*'         | 'hell*'
        'hell?-*'       | 'hell?-*'
        '\\*ello??'     | '\\Xello??'
        'hello\\[a-b\\]'| 'hello\\Xa-b\\X'
    }


}
