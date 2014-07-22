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

package nextflow.file.ggfs
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributeView
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.spi.FileSystemProvider

import nextflow.Session
import nextflow.extension.FilesExtensions
import nextflow.processor.TaskProcessor
import nextflow.script.BaseScript
import nextflow.file.FileHelper
import org.apache.commons.io.IOUtils
import org.gridgain.grid.GridGain
import org.gridgain.grid.logger.slf4j.GridSlf4jLogger
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GgFileSystemProviderTest extends Specification {

    @Shared
    def GgFileSystemProvider provider

    @Shared
    def GgFileSystem fs

    // run before the first feature method
    def setupSpec() {
        System.properties.setProperty('GRIDGAIN_HOME', new File('.').absolutePath)

        /*
         * configure and launch the grid
         */
        def file = new File('subprojects/nxf-gridgain/src/test/ggfs-test.xml')
        if( !file.exists() )
            file = new File('src/test/ggfs-test.xml')

        assert file.exists()
        def cfg = GridGain.loadConfiguration(file.absolutePath).get1()
        cfg.setGridLogger( new GridSlf4jLogger() )
        cfg.setGridGainHome( new File('.').absolutePath )
        def grid = GridGain.start(cfg)

        /*
         * get the installed provider for GgFileSystem
         */
        provider = (GgFileSystemProvider)FileSystemProvider.installedProviders().find { it.scheme == 'ggfs' }
        fs = (GgFileSystem)provider.newFileSystem( grid:grid )

    }


    def testGetFileSystem() {
        expect:
        provider.getFileSystem(URI.create('ggfs://x/y/z')) == fs
    }

    def testGetPath() {
        expect:
        provider.getPath(URI.create('ggfs:///path/to/file')) == new GgPath(fs, '/path/to/file')
    }

    def testGetScheme() {
        expect:
        provider.getScheme() == 'ggfs'
    }

    def testIsSameFile() {

        expect:
        provider.isSameFile(new GgPath(fs, '/file/x/y'), new GgPath(fs, '/file/x/y'))
        provider.isSameFile(new GgPath(fs, '/file/x/y'), new GgPath(fs, 'file/x/y'))
        provider.isSameFile(new GgPath(fs, '/file/x/y'), new GgPath(fs, '/file/../file/x/y'))
        !provider.isSameFile(new GgPath(fs, '/file/x/y'), new GgPath(fs, 'file/x/z'))

    }

    static rnd = new Random()

    def rndName( String prefix = 'file') {
        prefix + rnd.nextInt(10_000)
    }

    def testFileCreate() {

        when:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        def path2 = new GgPath(fs, rndName('/x/y/z/file') )
        Files.createFile(path1)

        then:
        Files.exists(path1)
        !Files.exists(path2)

        when:
        Files.createFile(path1)
        then:
        thrown(FileAlreadyExistsException)

    }

    def testFileExists() {

        when:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        def path2 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        then:
        Files.exists(path1)
        !Files.exists(path2)
        Files.readAllBytes(path1) == 'Hello'.bytes

    }

    def testFileCopy() {
        given:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        def path2 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        when:
        provider.copy(path1, path2)

        then:
        Files.exists(path1)
        Files.exists(path2)
        path2.text == 'Hello'
        !provider.isSameFile(path1, path2)

    }

    def testFileCopy2() {

        given:
        def path1 = Files.createTempFile('test',null)
        def path2 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        when:
        Files.copy(path1, path2)

        then:
        Files.exists(path1)
        Files.exists(path2)
        path2.text == 'Hello'
        !provider.isSameFile(path1, path2)

        cleanup:
        path1?.delete()

    }

    def testFileCopy3() {
        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('file1').text = 'Ciao'
        folder.resolve('file2').text = 'Hello'
        folder.resolve('sub_dir').mkdir()
        folder.resolve('sub_dir').resolve('file3').text = 'Hola'

        def target =  new GgPath(fs, rndName('/test-folder') )

        when:
        def result = FilesExtensions.copyTo(folder, target)
        then:
        result.resolve('file1').text == 'Ciao'
        result.resolve('file2').text == 'Hello'
        result.resolve('sub_dir').resolve('file3').text == 'Hola'

        cleanup:
        folder?.deleteDir()

    }


    def testFileMove() {
        given:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        def path2 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        when:
        provider.move(path1, path2)

        then:
        !Files.exists(path1)
        Files.exists(path2)
        path2.text == 'Hello'

    }


    def testFileReadAttrsFile() {

        given:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        when:
        def attrs = provider.readAttributes(path1, BasicFileAttributes)

        then:
        attrs.creationTime() == null
        attrs.lastModifiedTime().toMillis() > 0 && attrs.lastModifiedTime().toMillis() <= System.currentTimeMillis()
        attrs.size() == 5
        !attrs.directory
        attrs.regularFile
        !attrs.isSymbolicLink()

    }

    def testFileReadAttrsView() {

        given:
        def path1 = new GgPath(fs, rndName('/x/y/z/file') )
        path1.text = 'Hello'

        when:
        def view = provider.getFileAttributeView(path1, BasicFileAttributeView)
        def attrs = view.readAttributes()

        then:
        view instanceof GgFileAttributeView
        attrs.creationTime() == null
        attrs.lastModifiedTime().toMillis() > 0 && attrs.lastModifiedTime().toMillis() <= System.currentTimeMillis()
        attrs.size() == 5
        !attrs.directory
        attrs.regularFile
        !attrs.isSymbolicLink()

    }

    def testFileReadAttrsDir() {

        given:
        def path1 = new GgPath(fs, rndName('/x/y/z/dir') )
        Files.createDirectories(path1)

        when:
        def attrs = provider.readAttributes(path1, BasicFileAttributes)

        then:
        attrs.creationTime() == null
        attrs.lastModifiedTime().toMillis() > 0 && attrs.lastModifiedTime().toMillis() <= System.currentTimeMillis()
        attrs.size() == 0
        attrs.directory
        !attrs.regularFile
        !attrs.isSymbolicLink()

    }


    def testDirectoryStream() {

        given:
        def dir1 = new GgPath(fs, rndName('/stream/dir') )
        def dir2 = new GgPath(fs, rndName('/stream/dir') )
        def file1 = new GgPath(fs, rndName('/stream/file') )
        def file2 = new GgPath(fs, rndName('/stream/file') )
        def file3 = new GgPath(fs, rndName('/stream/file') )
        def file4 = dir2.resolve(rndName('file'))

        Files.createDirectories(dir1)
        Files.createDirectories(dir2)

        file1.text = 'file 1'
        file2.text = 'file 2'
        file3.text = 'file 3'
        file4.text = 'file 4'

        when:
        List<GgPath> list = provider.newDirectoryStream(new GgPath(fs, '/stream' ), { true } ).iterator().collect { it }

        then:
        list.size() == 5
        list.contains(dir1)
        list.contains(dir2)
        list.contains(file1)
        list.contains(file2)
        list.contains(file3)
        !list.contains(file4)

    }

    def testNewInputStream() {

        given:
        def file1 = new GgPath(fs, rndName('/x/file') )
        file1.text = 'Hello\nworld!'

        expect:
        IOUtils.readLines(provider.newInputStream(file1)) == ['Hello','world!']

        when:
        provider.newInputStream(file1, StandardOpenOption.WRITE)
        then:
        thrown(UnsupportedOperationException)

    }

    def testNewOutputStream() {

        when:
        def file1 = new GgPath(fs, rndName('/x/file') )
        provider.newOutputStream(file1, StandardOpenOption.READ)
        then:
        thrown(IllegalArgumentException)

        when:
        file1.text = 'Hello'
        // by default, a write over an existing channel overwrite it
        def output = provider.newOutputStream(file1)
        output.write( 'Ciao'.bytes )
        output.close()
        then:
        file1.text == 'Ciao'

        when:
        def file2 = new GgPath(fs, rndName('/x/file') )
        file2.text = 'Hello'
        // append to the existing file
        output = provider.newOutputStream(file2, StandardOpenOption.APPEND)
        output.write( ' world'.bytes )
        output.close()
        then:
        file2.text == 'Hello world'

    }

    def testDelete() {

        given:
        def file1 = new GgPath(fs, rndName('/x/file') )
        file1.text = 'some content'

        expect:
        Files.exists(file1)

        when:
        Files.delete(file1)
        then:
        !Files.exists(file1)

    }


    def testTempFolder() {

        when:
        def file1 = new GgPath(fs, rndName('/work') )
        def folder = FileHelper.createTempFolder(file1)
        then:
        folder.exists()

        when:
        def file2 = folder.resolve('just_a_file')
        file2.text = 'Hello'
        then:
        file2.text == 'Hello'

    }

    def testTempFolder2() {

        given:
        def session = new Session()
        session.workDir = new GgPath(fs, rndName('/work') )
        def binding = new Binding()
        binding.setVariable('__$session', session)

        def script = new BaseScript(binding) { Object run() { return null } }
        def proc = [:] as TaskProcessor
        proc.ownerScript = script

        when:
        def holder = proc.normalizeInputToFile('Hello', 'input.1')
        println ">> ${holder.storePath}"
        then:
        holder.sourceObj == 'Hello'
        holder.storePath.text == 'Hello'
    }

    def testGetFileAttributeView() {

        given:
        def path = new GgPath(fs, rndName('/x/y/z/file') )

        expect:
        provider.getFileAttributeView(path,BasicFileAttributeView) instanceof GgFileAttributeView
        provider.getFileAttributeView(path,GgFileAttributeView) instanceof GgFileAttributeView

    }


    // run after the last feature method
    def cleanupSpec() {
        GridGain.stop(true)
    }


}
