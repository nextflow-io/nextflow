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

package nextflow
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Paths

import nextflow.util.ArrayTuple
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowTest extends Specification {

    @Requires({ System.getenv('CI_GROOVY_VERSION') })
    def 'should match CI groovy version'() {
        expect:
        System.getenv('CI_GROOVY_VERSION') == GroovySystem.getVersion()
    }

    def testFile() {

        expect:
        Nextflow.file('file.log').toFile() == new File('file.log').canonicalFile
        Nextflow.file('relative/file.test').toFile() == new File( new File('.').canonicalFile, 'relative/file.test')
        Nextflow.file('/user/home/file.log').toFile() == new File('/user/home/file.log')

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

    def testFile3() {
        Exception e
        when:
        Nextflow.file(null)
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'Argument of `file` function cannot be null'

        when:
        Nextflow.file('')
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'Argument of `file` function cannot be empty'
    }

    def testFileWithWildcards() {

        setup:
        def folder = Files.createTempDirectory('nxftest') .toAbsolutePath()
        folder.resolve('hola1').text = 'abc'
        folder.resolve('helo2').text = 'abc'
        folder.resolve('hello3').text = 'abc'

        expect:
        Nextflow.file("$folder/ciao*") == []
        Nextflow.file("$folder/hel*").collect { it.name }.sort() == [ 'hello3','helo2' ]
        Nextflow.file("$folder/hol??").collect { it.name } == [ 'hola1' ]


        cleanup:
        folder?.delete()

    }

    def testFileWithWildcards2() {

        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        Files.createSymbolicLink(folder.resolve('link_file2.fa'), folder.resolve('file2.fa'))
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file_3.txt').text = 'file 3'
        folder.resolve('dir1').resolve('file_4.txt').text = 'file 4'
        folder.resolve('dir1').resolve('dir2').mkdir()

        when:
        def result = Nextflow.files(folder.resolve('file1.txt').toString())
        then:
        result.collect { it.name }.sort() == ['file1.txt']

        when:
        result = Nextflow.files("$folder/file*")
        then:
        result.collect { it.name }.sort() == ['file1.txt', 'file2.fa']

        when:
        result = Nextflow.files("$folder/dir1/dir2")
        then:
        result.collect { it.name }.sort()== ['dir2']

        when:
        result = Nextflow.files("$folder/**/file_*")
        then:
        result.collect { it.name }.sort() == ['file_3.txt', 'file_4.txt']

        when:
        result = Nextflow.files("$folder/*.fa")
        then:
        result.collect { it.name }.sort() == ['file2.fa', 'link_file2.fa']

        when:
        result = Nextflow.files("$folder/{file1.txt,dir1/file_4.txt}")
        then:
        result.collect { it.name }.sort() == ['file1.txt','file_4.txt']

        when:
        result = Nextflow.files("$folder/*")
        then:
        result.collect { it.name }.sort() == ['file1.txt', 'file2.fa', 'link_file2.fa']

        when:
        result = Nextflow.files("$folder/*", type:'file')
        then:
        result.collect { it.name }.sort() == ['file1.txt', 'file2.fa', 'link_file2.fa']

        when:
        result = Nextflow.files("$folder/*", type:'dir')
        then:
        result.collect { it.name }.sort() == ['dir1']

        cleanup:
        folder?.deleteDir()

    }

    def testFileWithWildcards3() {

        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        Files.createSymbolicLink( folder.resolve('file_link.fa'), folder.resolve('file2.fa'))
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file3.txt').text = 'file 3'
        folder.resolve('dir1').resolve('dir2').mkdirs()
        folder.resolve('dir1').resolve('dir2').resolve('file4.fa').text = 'file '
        Files.createSymbolicLink( folder.resolve('dir_link'), folder.resolve('dir1') )

        when:
        def result = Nextflow.files("$folder/**.fa", relative: true)
        then:
        result.collect { it.toString() } .sort() == ['dir1/dir2/file4.fa', 'dir_link/dir2/file4.fa', 'file2.fa', 'file_link.fa']

        when:
        result = Nextflow.files("$folder/**.fa", relative: true, followLinks: false)
        then:
        result.collect { it.toString() } .sort() == ['dir1/dir2/file4.fa', 'file2.fa']

        when:
        result = Nextflow.files("$folder/**.fa", relative: true, maxDepth: 1)
        then:
        result.collect { it.toString() } .sort() == ['file2.fa', 'file_link.fa']

        when:
        result = Nextflow.files("$folder/**", relative: true, type:'file')
        then:
        result.collect { it.toString() } .sort() == ['dir1/dir2/file4.fa',
                                                     'dir1/file3.txt',
                                                     'dir_link/dir2/file4.fa',
                                                     'dir_link/file3.txt',
                                                     'file1.txt',
                                                     'file2.fa',
                                                     'file_link.fa']

        when:
        result = Nextflow.files("$folder/**", relative: true, type:'file', followLinks: false)
        then:
        result.collect { it.toString() } .sort() == ['dir1/dir2/file4.fa',
                                                     'dir1/file3.txt',
                                                     'file1.txt',
                                                     'file2.fa']

        when:
        result = Nextflow.files("$folder/**", relative: true, type:'dir')
        then:
        result.collect { it.toString() } .sort() == ['dir1',
                                                     'dir1/dir2',
                                                     'dir_link',
                                                     'dir_link/dir2' ]


        when:
        result = Nextflow.files("$folder/**", relative: true, type:'dir', followLinks: false)
        then:
        result.collect { it.toString() } .sort() == ['dir1', 'dir1/dir2']


        when:
        result = Nextflow.files("$folder/**", relative: true, type:'any', followLinks: false)
        then:
        result.collect { it.toString() } .sort() == ['dir1',
                                                     'dir1/dir2',
                                                     'dir1/dir2/file4.fa',
                                                     'dir1/file3.txt',
                                                     'file1.txt',
                                                     'file2.fa']

        cleanup:
        folder?.deleteDir()
    }

    def testVisitNoSuchFile() {

        when:
        def result = Nextflow.file('/some/missing/path/*')
        then:
        result == []
    }


    def testTuple() {

        expect:
        Nextflow.tuple(1) == new ArrayTuple([1])
        Nextflow.tuple(1,2,3) == new ArrayTuple([1,2,3])
        Nextflow.tuple([4,5,6]) == new ArrayTuple([4,5,6])
        Nextflow.tuple([4,5,6], [1,2]) == new ArrayTuple([[4,5,6], [1,2]])
    }


    def 'should match file with glob pattern' () {

        given:
        def root = Files.createTempDirectory('test')
        def file1 = root.resolve('file-*.txt')
        def file2 = root.resolve('file-?.txt')
        def file3 = root.resolve('file[a-b].txt')
        def file4 = root.resolve('file{a,b}.txt')
        def filea = root.resolve('filea.txt')
        def fileb = root.resolve('fileb.txt')

        when:
        file3.text = 'Hello'
        file4.text = 'World'
        file1.text = 'Star'
        file2.text = 'Mark'
        filea.text = 'aaa'
        fileb.text = 'bbb'

        then:
        file1.exists()
        file2.exists()
        file3.exists()
        file4.exists()
        filea.exists()
        fileb.exists()
        file1.name == 'file-*.txt'
        file2.name == 'file-?.txt'
        file3.name == 'file[a-b].txt'
        file4.name == 'file{a,b}.txt'

        expect:
        Nextflow.file(root.resolve('file-\\*.txt')) == file1
        Nextflow.file(root.resolve('file-\\?.txt')) == file2
        Nextflow.file(root.resolve('file-*')) *. name .sort() == ['file-*.txt', 'file-?.txt']
        Nextflow.file(root.resolve('file[a-b].*')) *. name .sort() == ['filea.txt', 'fileb.txt']

        Nextflow.file(root.resolve('file-*.txt'), glob: false) == file1
        Nextflow.file("$root/file-*.txt", glob: false) == file1

        cleanup:
        root?.deleteDir()
    }

    def 'should check if file exists' () {
        given:
        def folder = Files.createTempDirectory('test')
        def foo = folder.resolve('foo.txt')

        when:
        Nextflow.file(foo.toString(), checkIfExists: true)
        then:
        def e = thrown(NoSuchFileException)
        e.message == foo.toString()

        cleanup:
        folder?.deleteDir()
    }

}
