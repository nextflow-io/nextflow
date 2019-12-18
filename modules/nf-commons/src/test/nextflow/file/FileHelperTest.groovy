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

package nextflow.file

import spock.lang.Unroll

import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.nio.file.spi.FileSystemProvider

import com.google.common.jimfs.Configuration
import com.google.common.jimfs.Jimfs
import nextflow.Global
import nextflow.ISession
import spock.lang.Specification

import static java.nio.file.LinkOption.NOFOLLOW_LINKS
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperTest extends Specification {

    static private fs = Jimfs.newFileSystem(Configuration.unix());

    static Path createInMemTempDir() {
        Path tmp = fs.getPath("/tmp");
        tmp.mkdir()
        Files.createTempDirectory(tmp, 'test')
    }

    def 'should return a Path object' () {

        expect:
        FileHelper.asPath('file.txt') == Paths.get('file.txt')
        FileHelper.asPath('file:file.txt') == Paths.get('file:file.txt')
        FileHelper.asPath('file:/some/file.txt') == Paths.get('/some/file.txt')
        FileHelper.asPath('file:///some/file.txt') == Paths.get('/some/file.txt')
        FileHelper.asPath('file:///some/fil?.{fa,txt}') == Paths.get('/some/fil?.{fa,txt}')

        when:
        FileHelper.asPath('file://some/file.txt')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should create a valid uri' () {

        def URI uri

        when:
        uri = FileHelper.toPathURI('jimfs:test/some/file')
        then:
        uri == new URI('jimfs:test/some/file')
        uri.scheme == 'jimfs'
        uri.path == null        // note: path is null because it is not a valid hierarchical URI. See http://docs.oracle.com/javase/7/docs/api/java/net/URI.html

        when:
        uri = FileHelper.toPathURI('/a/b/c')
        then:
        uri == new URI('/a/b/c')
        uri.path == '/a/b/c'
        uri.authority == null

        when:
        uri = FileHelper.toPathURI('a/b/c')
        then:
        uri == new URI('a/b/c')
        uri.path == 'a/b/c'
        uri.authority == null

        when:
        uri = FileHelper.toPathURI('//xxx/a/b/c')
        then:
        uri == new URI('//xxx/a/b/c')
        uri.authority == 'xxx'
        uri.path == '/a/b/c'

        when:
        uri = FileHelper.toPathURI('file:/a/b/c')
        then:
        uri == new URI('file:/a/b/c')
        uri.path == '/a/b/c'
        uri.scheme == 'file'

        when:
        uri = FileHelper.toPathURI('file://a/b/c')
        then:
        uri == new URI('file://a/b/c')
        uri.path == '/b/c'
        uri.scheme == 'file'
        uri.authority == 'a'

        when:
        uri = FileHelper.toPathURI('file:///a/b/c')
        then:
        uri == new URI('file:///a/b/c')
        uri.path == '/a/b/c'
        uri.scheme == 'file'

        when:
        uri = FileHelper.toPathURI('file:///a/b/c{1,2}')
        then:
        uri.path == '/a/b/c{1,2}'
        uri.scheme == 'file'

        when:
        uri = FileHelper.toPathURI('s3:///cbcrg-eu/raw/**_R1*{fastq,fq,fastq.gz,fq.gz}')
        then:
        uri.path == '/cbcrg-eu/raw/**_R1*{fastq,fq,fastq.gz,fq.gz}'
        uri.scheme == 's3'

        when:
        uri = FileHelper.toPathURI('s3:///cbcrg-eu//raw/x_r1.fq')
        then:
        uri == new URI('s3:///cbcrg-eu//raw/x_r1.fq')
        uri.path == '/cbcrg-eu//raw/x_r1.fq'
        uri.scheme == 's3'

        when:
        uri = FileHelper.toPathURI('s3://cbcrg-eu//raw/x_r1.fq')
        then:
        uri == new URI('s3:///cbcrg-eu//raw/x_r1.fq')
        uri.path == '/cbcrg-eu//raw/x_r1.fq'
        uri.scheme == 's3'

        when:
        uri = FileHelper.toPathURI('dxfs://grape:/data/ggal/ggal_test_1.fq')
        then:
        uri == new URI('dxfs://grape:/data/ggal/ggal_test_1.fq')
        uri.scheme == 'dxfs'
        uri.path == '/data/ggal/ggal_test_1.fq'
        uri.authority == 'grape:'

        when:
        uri = FileHelper.toPathURI('dxfs://grape:/data/ggal/ggal_test_?.fq')
        then:
        uri == new URI('dxfs://grape:/data/ggal/ggal_test_%3F.fq')
        uri.scheme == 'dxfs'
        uri.path == '/data/ggal/ggal_test_?.fq'
        uri.authority == 'grape:'

        when:
        uri = FileHelper.toPathURI('file:///some/file.txt#abc')
        then:
        uri.scheme == 'file'
        uri.path == '/some/file.txt#abc'

        when:
        uri = FileHelper.toPathURI('igfs:///work/dir')
        then:
        uri.scheme == 'igfs'
        uri.path == '/work/dir'

        when:
        uri = FileHelper.toPathURI('igfs://work/dir')
        then:
        uri.scheme == 'igfs'
        uri.path == '/work/dir'

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


    def 'get env map'() {

        given:
        def sess = Global.session = Mock(ISession)
        def env = [:]
        env.put('AWS_ACCESS_KEY','a1')
        env.put('AWS_SECRET_KEY','s1')

        expect:
        // properties have priority over the environment map
        FileHelper.envFor0('s3', env).access_key == 'a1'
        FileHelper.envFor0('s3', env).secret_key == 's1'

        // any other return just the session
        FileHelper.envFor0('dxfs', env).session == sess

    }

    def 'cached path'() {

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


    def 'visit files'() {

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

    def 'visit files 2'() {

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

    def 'visit files 3' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'
        folder.resolve('file3.bam').text = 'file 2'
        Files.createSymbolicLink(
                folder.resolve('file4.fa'),
                folder.resolve('file3.bam'))

        when:
        def result = []
        FileHelper.visitFiles(folder, '*.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file2.fa','file4.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, '*.fa', relative: true, followLinks: false) { result << it.toString() }
        then:
        result.sort() == ['file2.fa']

        cleanup:
        folder?.deleteDir()
    }

    def 'visit files in a base path with glob characters' () {
        given:
        def folder = Files.createTempDirectory('test[a-b]')
        folder.resolve('file1.txt').text = 'file 1'
        folder.resolve('file2.fa').text = 'file 2'

        when:
        def result = []
        FileHelper.visitFiles(folder, 'file*', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file1.txt', 'file2.fa']

        cleanup:
        folder.deleteDir()
    }

    def 'match files containing glob characters' () {
        given:
        def folder = Files.createTempDirectory('test[a-b]')
        folder.resolve('file*.txt').text = 'file 1'
        folder.resolve('file?.fa').text = 'file 2'
        folder.resolve('file{a,b}.fa').text = 'file 3'
        folder.resolve('file[a-b].fa').text = 'file 4'

        when:
        def result = []
        FileHelper.visitFiles(folder, 'x*', relative: true) { result << it.toString() }
        then:
        result.sort() == []

        when:
        result = []
        FileHelper.visitFiles(folder, 'file\\*.txt', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file*.txt']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file\\?.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file?.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file\\{a,b\\}.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file{a,b}.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file\\[a-b\\].fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file[a-b].fa']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file*.fa', relative: true) { result << it.toString() }
        then:
        result.sort() == ['file?.fa','file[a-b].fa','file{a,b}.fa']

        when:
        result = []
        FileHelper.visitFiles(folder, 'file\\[*\\].fa', relative: true) { result << it.toString() }
        then:
        result.sort() == [ 'file[a-b].fa' ]


        cleanup:
        folder.deleteDir()
    }

    def 'get max depth'() {
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


    def 'no such file'() {

        when:
        FileHelper.visitFiles(Paths.get('/some/missing/path'),'*', { return it })
        then:
        thrown(NoSuchFileException)

    }


    def 'should copy file path to foreign file system' () {

        given:
        def target = createInMemTempDir()
        def source = Files.createTempDirectory('test')
        def file = source.resolve('file.txt')
        def copy = target.resolve('copy.txt')
        file.text = 'Hello world'

        when:
        FileHelper.copyPath(file, copy)
        then:
        copy.text == 'Hello world'
        file.exists()

        when:
        FileHelper.copyPath(file, copy)
        then:
        thrown(FileAlreadyExistsException)

        when:
        file.text = 'Hi there'
        FileHelper.copyPath(file, copy, StandardCopyOption.REPLACE_EXISTING)
        then:
        copy.text == 'Hi there'
        file.exists()

        cleanup:
        source.deleteDir()

    }

    def 'should copy dir path to foreign file system' () {

        given:
        def root = Files.createTempDirectory('test')
        def source = root.resolve('source')
        def target = createInMemTempDir().resolve('copy')

        new FileTreeBuilder(root.toFile())
                .dir('source') {
                    file('file_1','alpha')
                    file('file_2','beta')
                    dir('dir1')  {
                        file('file_3', 'delta')
                    }
                    dir('dir2') {
                        file('file_4', 'gamma')
                        dir('dir3') {
                            file('file_5', 'omega')
                        }
                    }
                }

        when:
        FileHelper.copyPath(source, target)
        then:
        target.isDirectory()
        target.resolve('file_1').text == 'alpha'
        target.resolve('file_2').text == 'beta'
        target.resolve('dir1/file_3').text == 'delta'
        target.resolve('dir2/file_4').text == 'gamma'
        target.resolve('dir2/dir3/file_5').text == 'omega'
        source.exists()

        when:
        target = root.resolve('copy')
        FileHelper.copyPath(source, target)
        then:
        target.isDirectory()
        target.resolve('file_1').text == 'alpha'
        target.resolve('file_2').text == 'beta'
        target.resolve('dir1/file_3').text == 'delta'
        target.resolve('dir2/file_4').text == 'gamma'
        target.resolve('dir2/dir3/file_5').text == 'omega'
        source.exists()

        root?.deleteDir()
    }

    def 'should copy path content' () {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('source').mkdir()
        folder.resolve('target').mkdir()

        folder.resolve('source/file.txt').text = 'hello world'
        folder.resolve('source/file.txt').mklink(folder.resolve('source/link.txt'))
        folder.resolve('source/dir-x').mkdir()
        folder.resolve('source/dir-x/content.txt').text = 'foo'
        folder.resolve('source/dir-x').mklink(folder.resolve('source/link-dir'))

        when:
        def s = folder.resolve('source/link.txt')
        def t = folder.resolve('target/file-b.txt')
        FileHelper.copyPath(s,t)
        then:
        t.text == 'hello world'
        Files.isRegularFile(t)

        when:
        s = folder.resolve('source/link.txt')
        t = folder.resolve('target/file-c.txt')
        FileHelper.copyPath(s,t, LinkOption.NOFOLLOW_LINKS)
        then:
        t.text == 'hello world'
        Files.isSymbolicLink(t)

        when:
        s = folder.resolve('source/link-dir')
        t = folder.resolve('target/copy-dir')
        FileHelper.copyPath(s,t)
        then:
        Files.isDirectory(t)
        t.resolve('content.txt').text == 'foo'

        when:
        s = folder.resolve('source/link-dir')
        t = folder.resolve('target/link-dir')
        FileHelper.copyPath(s,t, LinkOption.NOFOLLOW_LINKS)
        then:
        Files.isSymbolicLink(t)
        Files.isDirectory(t)
        t.resolve('content.txt').text == 'foo'

        cleanup:
        folder?.deleteDir()
    }

    def 'should move file path to foreign file system' () {

        given:
        def target = createInMemTempDir()
        def source = Files.createTempDirectory('test')
        def file = source.resolve('file.txt')
        def copy = target.resolve('copy.txt')
        file.text = 'Hello world'

        when:
        FileHelper.movePath(file, copy)
        then:
        copy.text == 'Hello world'
        !file.exists()

    }

    def 'should move dir path to foreign file system' () {

        given:
        def root = Files.createTempDirectory('test')
        def source = root.resolve('source')
        def target = createInMemTempDir().resolve('copy')

        new FileTreeBuilder(root.toFile())
                .dir('source') {
                    file('file_1','alpha')
                    file('file_2','beta')
                    dir('dir1')  {
                        file('file_3', 'delta')
                    }
                    dir('dir2') {
                        file('file_4', 'gamma')
                        dir('dir3') {
                            file('file_5', 'omega')
                        }
                    }
        }

        when:
        FileHelper.movePath(source, target)
        then:
        target.isDirectory()
        target.resolve('file_1').text == 'alpha'
        target.resolve('file_2').text == 'beta'
        target.resolve('dir1/file_3').text == 'delta'
        target.resolve('dir2/file_4').text == 'gamma'
        target.resolve('dir2/dir3/file_5').text == 'omega'
        !source.exists()

        cleanup:
        root?.deleteDir()

    }

    def 'should move dir path to local file system' () {

        given:
        def root = Files.createTempDirectory('test')
        def source = root.resolve('source')
        def target = root.resolve('copy')

        new FileTreeBuilder(root.toFile())
                .dir('source') {
                    file('file_1','alpha')
                    file('file_2','beta')
                    dir('dir1')  {
                        file('file_3', 'delta')
                    }
                    dir('dir2') {
                        file('file_4', 'gamma')
                        dir('dir3') {
                            file('file_5', 'omega')
                        }
                    }
        }

        when:
        FileHelper.movePath(source, target)
        then:
        target.isDirectory()
        target.resolve('file_1').text == 'alpha'
        target.resolve('file_2').text == 'beta'
        target.resolve('dir1/file_3').text == 'delta'
        target.resolve('dir2/file_4').text == 'gamma'
        target.resolve('dir2/dir3/file_5').text == 'omega'
        !source.exists()

        cleanup:
        root?.deleteDir()

    }

    def 'should delete file or a directory'() {

        given:
        def folder = Files.createTempDirectory('test')

        def file_1 = folder.resolve('file_1')
        def file_2 = folder.resolve('file_2')
        def dir_1 = folder.resolve('dir_1')
        def dir_2 = folder.resolve('dir_2')
        def link_to_dir_2 = folder.resolve('link_to_dir_2')
        def link_to_file_2 = folder.resolve('link_to_file_2')
        def link_to_missing = folder.resolve('link_to_missing')

        Files.createFile(file_1)
        Files.createFile(file_2)
        Files.createDirectory(dir_1)
        Files.createDirectory(dir_2)
        Files.createSymbolicLink(link_to_file_2, file_2)
        Files.createSymbolicLink(link_to_dir_2, dir_2)
        Files.createSymbolicLink(link_to_missing, folder.resolve('missing'))

        Files.createFile(dir_1.resolve('hello'))
        Files.createDirectory(dir_1.resolve('sub_dir'))
        Files.createFile(dir_1.resolve('sub_dir/world'))

        expect:
        FileHelper.deletePath(file_1)
        // delete the file
        !Files.exists(file_1)

        FileHelper.deletePath(dir_1)
        // delete directory including sub directory
        !Files.exists(dir_1)

        FileHelper.deletePath(link_to_file_2)
        // delete the symlink but not the target file
        !Files.exists(link_to_file_2 )
        !Files.exists(link_to_file_2, LinkOption.NOFOLLOW_LINKS)
        Files.exists(file_2)

        FileHelper.deletePath(link_to_dir_2)
        // delete the symlink but not the target dir
        !Files.exists(link_to_dir_2)
        !Files.exists(link_to_dir_2, LinkOption.NOFOLLOW_LINKS)
        Files.exists(dir_2)     // <-- note: delete the symlink but not the target file

        FileHelper.deletePath(link_to_missing)
        // delete the symlink
        !Files.exists(link_to_missing)
        !Files.exists(link_to_missing, LinkOption.NOFOLLOW_LINKS)

        cleanup:
        folder?.deleteDir()

    }

    def 'should read file attributes' () {
        given:
        def folder = Files.createTempDirectory('test')

        def file_1 = folder.resolve('file_1.txt')
        def dir_1 = folder.resolve('dir_1')
        def link_to_file = folder.resolve('link_1')
        def link_to_dir = folder.resolve('link_2')

        when:
        Files.createFile(file_1)
        then:
        FileHelper.readAttributes(file_1).isRegularFile()

        when:
        Files.createSymbolicLink(link_to_file, file_1)
        then:
        FileHelper.readAttributes(link_to_file).isRegularFile()
        FileHelper.readAttributes(link_to_file,NOFOLLOW_LINKS).isSymbolicLink()

        when:
        Files.createDirectory(dir_1)
        then:
        FileHelper.readAttributes(dir_1).isDirectory()

        when:
        Files.createSymbolicLink(link_to_dir, dir_1)
        then:
        FileHelper.readAttributes(link_to_dir).isDirectory()
        FileHelper.readAttributes(link_to_dir,NOFOLLOW_LINKS).isSymbolicLink()

        when:
        Files.delete(file_1)
        then:
        FileHelper.readAttributes(link_to_file) == null

        cleanup:
        folder?.deleteDir()
    }

    def 'should check if glob is allowed' () {

        given:
        def FS = Mock(FileSystem)
        def PROVIDER = Mock(FileSystemProvider)
        def PATH = Mock(Path)

        when:
        def result = FileHelper.isGlobAllowed(PATH)
        then:
        PATH.getFileSystem() >> FS
        FS.provider() >> PROVIDER
        PROVIDER.getScheme() >> SCHEME
        result == EXPECTED
        
        where:
        SCHEME      | EXPECTED
         'file'     |  true
         's3'       |  true
         'http'     |  false
         'https'    |  false
         'ftp'      |  false
         'ftps'     |  false
    }

    @Unroll
    def 'should return an identifier given file #PATH'() {

        expect:
        FileHelper.getIdentifier(Paths.get(PATH)) == EXPECTED

        where:
        PATH                | EXPECTED
        'foo/bar.nf.txt'    | 'bar'
        'foo-bar-baz.nf'    | 'foo_bar_baz'
        '123-fo0'           | '_23_fo0'
        '--a  b  c'         | '_a_b_c'
    }

    @Unroll
    def 'should get url protocol' () {
        expect:
        FileHelper.getUrlProtocol(STR)  == EXPECTED
        where:
        EXPECTED    | STR
        'ftp'       | 'ftp://abc.com'
        's3'        | 's3://bucket/abc'
        null        | '3s://bucket/abc'
        null        | 'abc:xyz'
        null        | '/a/bc/'
    }

}
