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

package nextflow.file
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths

import com.google.common.jimfs.Configuration
import com.google.common.jimfs.Jimfs
import nextflow.Global
import nextflow.ISession
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperTest extends Specification {

    def 'test asPath' () {

        given:
        Jimfs.newFileSystem(Configuration.unix());

        expect:
        FileHelper.asPath('file.txt') == Paths.get('file.txt')
        FileHelper.asPath('file:///file.txt') == Paths.get( URI.create('file:///file.txt') )
        FileHelper.asPath('jimfs:test/some/file') == Paths.get('jimfs:test/some/file')

    }

    def 'test normalize' () {

        expect:
        FileHelper.normalizePath( 'file.name' ) == 'file.name'
        FileHelper.normalizePath( '~file.name' ) == '~file.name'
        FileHelper.normalizePath( '~' ) == System.properties['user.home']
        FileHelper.normalizePath( '~/path/file.name' ) == System.properties['user.home'] + '/path/file.name'

    }

    def 'test isEmpty file'() {

        setup:
        def emptyFile = File.createTempFile('test','test')
        def notEmptyFile = File.createTempFile('test','test')
        notEmptyFile.text = 'HOLA'
        emptyFile.deleteOnExit()
        notEmptyFile.deleteOnExit()

        expect:
        FileHelper.empty(emptyFile)
        !FileHelper.empty(notEmptyFile)

    }

    def 'test isEmpty dir'() {

        setup:
        def emptyDir = File.createTempDir()
        def notEmptyDir = File.createTempDir()
        File.createTempFile('test','test', notEmptyDir)

        expect:
        FileHelper.empty(emptyDir)
        !FileHelper.empty(notEmptyDir)

        cleanup:
        emptyDir.deleteDir()
        notEmptyDir.deleteDir()

    }


    def 'test empty file' () {

        setup:
        def fileEmpty = File.createTempFile('test','test')
        fileEmpty.deleteOnExit()

        File folderEmpty = File.createTempDir()

        File  folderNotEmpty = File.createTempDir()
        def fileInFolder = new File(folderNotEmpty, 'filename')
        fileInFolder.createNewFile()

        def fileNotEmpty = File.createTempFile('test','test')
        fileNotEmpty.text = 'Hola'
        fileNotEmpty.deleteOnExit()

        expect:
        FileHelper.empty(new File('non existing'))
        FileHelper.empty(fileEmpty)
        !FileHelper.empty(fileNotEmpty)
        FileHelper.empty(folderEmpty)
        !FileHelper.empty(folderNotEmpty)

        cleanup:
        fileEmpty.delete()
        folderNotEmpty?.deleteDir()
        folderEmpty?.deleteDir()

    }

    def 'test empty path' () {

        given:
        Path baseFolder = Files.createTempDirectory('empty')

        Path fileEmpty = Files.createTempFile(baseFolder, 'test','txt')
        Path folderEmpty = Files.createTempDirectory(baseFolder, null)
        Path folderNotEmpty = Files.createTempDirectory(baseFolder, null)
        Path fileInFolder = folderNotEmpty.resolve( 'empty_file' )
        Files.createFile(fileInFolder)

        Path fileNotEmpty = folderNotEmpty.resolve('not_empty_file')
        fileNotEmpty.text = 'Hola'

        Path fileNotExist = Paths.get('not_existing_file')

        expect:
        FileHelper.empty(fileNotExist)
        FileHelper.empty(fileEmpty)
        !FileHelper.empty(fileNotEmpty)
        FileHelper.empty(folderEmpty)
        !FileHelper.empty(folderNotEmpty)

        cleanup:
        baseFolder?.deleteDir()

    }



    def 'test nameParts' () {

        expect:
        FileHelper.nameParts("hola") == ['hola',0]
        FileHelper.nameParts("hola123")  == ['hola',123]
        FileHelper.nameParts("hola1")  == ['hola',1]
        FileHelper.nameParts("x")  == ['x',0 ]
        FileHelper.nameParts("x1")  == ['x',1 ]
        FileHelper.nameParts("1")  == ['',1]

    }



    def 'tests parent' () {

        expect:
        new File(".") .parentFile == null
        new File("hola") .parentFile == null
        new File("/").parentFile == null
        new File("/").absoluteFile == new File("/")
        new File('.').absolutePath.startsWith('/')
    }


    def testGeEnvMap() {

        given:
        def sess = Mock(ISession)
        Global.session = sess
        def env = [:]
        env.put('AWS_ACCESS_KEY','a1')
        env.put('AWS_SECRET_KEY','s1')

        expect:
        // properties have priority over the environment map
        filehelper.getenvmap0('s3', env) == [access_key:'a1', secret_key:'s1']
        // none of them
        FileHelper.getEnvMap0('s3', [:]) == [:]

        // any other return just the session
        FileHelper.getEnvMap0('dxfs', env).session == sess

    }

    def testCachedPath() {

        given:
        final localCacheDir = Files.createTempDirectory('cache')
        final testFolder = Files.createTempDirectory('test')
        final file1 = testFolder.resolve('file1.txt')
        final file2 = testFolder.resolve('file2.txt')
        file1.text = 'File 1'
        file2.text = 'File 2'

        // create a sub-dir with two files
        final dir1 = testFolder.resolve('dir1')
        dir1.mkdir()
        final file3 = dir1.resolve('file3')
        final file4 = dir1.resolve('file4')
        file3.text = 'File 3'
        file4.text = 'File 4'

        final id = UUID.randomUUID()

        when:
        def cache1 = FileHelper.getLocalCachePath(file1, localCacheDir, id)
        then:
        cache1 != file1
        cache1.text == file1.text
        cache1.getFileName() == file1.getFileName()
        cache1.parent.parent.parent == localCacheDir

        when:
        def cache2 = FileHelper.getLocalCachePath(file2, localCacheDir, id)
        then:
        cache2 != cache1
        cache2 != file2
        cache2.text == file2.text
        cache2.getFileName() == file2.getFileName()
        cache2.parent.parent.parent == localCacheDir

        when:
        def cache3 = FileHelper.getLocalCachePath(file1, localCacheDir, id)
        then:
        cache3 == cache1
        cache3 != file1
        cache3.text == file1.text

        when:
        def cache4 = FileHelper.getLocalCachePath(dir1, localCacheDir, id)
        then:
        Files.isDirectory(cache4)
        cache4 != dir1
        cache4.getName() == dir1.getName()
        cache4.resolve('file3').text == file3.text
        cache4.resolve('file4').text == file4.text

        cleanup:
        testFolder?.deleteDir()
        localCacheDir?.deleteDir()
    }


    def testVisitFiles() {

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
        def result = []
        FileHelper.visitFiles(folder, 'file1.txt', relative: true) { result << it.toString() }
        then:
        result == ['file1.txt']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file*', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file1.txt', 'file2.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, 'dir1/dir2', relative: true) { result << it.toString() }
        then:
        result.sort() == ['dir1/dir2']

        when:
        result = []
        FileHelper.visitFiles(folder, '**/file_*', relative: true) { result << it.toString() }
        then:
        result.sort() == ['dir1/file_3.txt', 'dir1/file_4.txt']

        when:
        result = []
        FileHelper.visitFiles(folder, '*.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file2.fa', 'link_file2.fa']

        cleanup:
        folder?.deleteDir()

    }

    def testVisitFiles2() {

        given:
        def folder = Files.createTempDirectory('test')

        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        folder.resolve('dir1').mkdir()
        folder.resolve('dir1').resolve('file3.txt').text = 'file 3'
        folder.resolve('dir1').resolve('dir2').mkdirs()
        folder.resolve('dir1').resolve('dir2').resolve('file4.fa').text = 'file '
        Files.createSymbolicLink( folder.resolve('dir_link'), folder.resolve('dir1') )

        when:
        def result = []
        FileHelper.visitFiles(folder, '**.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['dir1/dir2/file4.fa', 'dir_link/dir2/file4.fa', 'file2.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, '**.fa', relative: true, followLinks: false) { result << it.toString() }
        then:
        result.sort() == ['dir1/dir2/file4.fa', 'file2.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, '**.fa', relative: true, maxDepth: 1) { result << it.toString() }
        then:
        result.sort() == ['file2.fa']

        cleanup:
        folder?.deleteDir()
    }

    def testGetMAxDepth() {
        expect:
        FileHelper.getMaxDepth(1,null) == 1
        FileHelper.getMaxDepth(10,null) == 10
        FileHelper.getMaxDepth(0,'**') == 0
        FileHelper.getMaxDepth(1,'**') == 1
        FileHelper.getMaxDepth(null,null) == 0
        FileHelper.getMaxDepth(null,'abc') == 0
        FileHelper.getMaxDepth(null,'a/b') == 1
        FileHelper.getMaxDepth(null,'a/b/c') == 2
        FileHelper.getMaxDepth(null,'a/b/c/d') == 3
        FileHelper.getMaxDepth(null,'/a/b/c/d') == 3
        FileHelper.getMaxDepth(null,'a/**') == Integer.MAX_VALUE
    }

    def testGetPathAndPattern () {

        expect:
        FileHelper.getFolderAndPattern( '/some/file/name.txt' ) == ['/some/file/', 'name.txt', null]
        FileHelper.getFolderAndPattern( '/some/file/na*.txt' ) == ['/some/file/', 'na*.txt', null]
        FileHelper.getFolderAndPattern( '/some/file/na??.txt' ) == ['/some/file/', 'na??.txt', null]
        FileHelper.getFolderAndPattern( '/some/file/*.txt' ) == ['/some/file/', '*.txt', null]
        FileHelper.getFolderAndPattern( '/some/file/?.txt' ) == ['/some/file/', '?.txt', null]
        FileHelper.getFolderAndPattern( '/some/file/*' ) == ['/some/file/', '*', null]
        FileHelper.getFolderAndPattern( '/some/file/' ) == ['/some/file/', '', null]
        FileHelper.getFolderAndPattern( 'path/filename.txt' ) == ['path/', 'filename.txt', null]
        FileHelper.getFolderAndPattern( 'filename.txt' ) == ['./', 'filename.txt', null]
        FileHelper.getFolderAndPattern( './file.txt' ) == ['./', 'file.txt', null]

        FileHelper.getFolderAndPattern( '/some/file/**/*.txt' ) == ['/some/file/', '**/*.txt', null]

        FileHelper.getFolderAndPattern( 'dxfs:///some/file/**/*.txt' ) == ['/some/file/', '**/*.txt', 'dxfs']
        FileHelper.getFolderAndPattern( 'dxfs://some/file/**/*.txt' ) == ['some/file/', '**/*.txt', 'dxfs']
        FileHelper.getFolderAndPattern( 'dxfs://*.txt' ) == ['./', '*.txt', 'dxfs']
        FileHelper.getFolderAndPattern( 'dxfs:///*.txt' ) == ['/', '*.txt', 'dxfs']
        FileHelper.getFolderAndPattern( 'dxfs:///**/*.txt' ) == ['/', '**/*.txt', 'dxfs']

        FileHelper.getFolderAndPattern( 'file{a,b}') == ['./', 'file{a,b}', null]
        FileHelper.getFolderAndPattern( 'test/data/file{a,b}') == ['test/data/', 'file{a,b}', null]
        FileHelper.getFolderAndPattern( 'test/{file1,file2}') == ['test/', '{file1,file2}', null]
        FileHelper.getFolderAndPattern( '{file1,file2}') == ['./', '{file1,file2}', null]
        FileHelper.getFolderAndPattern( '{test/file1,data/file2}') == ['./', '{test/file1,data/file2}', null]
        FileHelper.getFolderAndPattern( 'data/{p/file1,q/file2}') == ['data/', '{p/file1,q/file2}', null]
    }

    def testNoSuchFile() {

        when:
        FileHelper.visitFiles(Paths.get('/some/missing/path'),'*', { return it })
        then:
        thrown(NoSuchFileException)

    }

    def 'test isGlobPattern' () {

        expect:
        FileHelper.isGlobPattern(pattern) == result

        where:
        pattern     | result
        'hola'      | false
        '1-2-3'     | false
        'hello.txt' | false
        'hello{x'   | false
        'hello[x'   | false
        'hello(a)'  | false
        'some/path' | false
        '*'         | true
        'hola*'     | true
        'hola?'     | true
        '?'         | true
        'hola[a]'   | true
        'hola[a-z]' | true
        'hola{a,b}' | true
        'hola{}'    | false
        'hola{a}'   | false
        'hola[]'    | false

    }

}
