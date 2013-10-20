package nextflow.extension

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FilesExtensionsTest extends Specification {



    def testGetExtension() {

        expect:
        new File('file.txt').getExtension() == 'txt'
        Paths.get('file.txt').getExtension() == 'txt'
        new File('file_name').getExtension() == ''
        Paths.get('file_name').getExtension() == ''

    }

    def testBaseName() {

        expect:
        new File('/path/file.txt').getBaseName() == 'file'
        Paths.get('/path/file.txt').getBaseName() == 'file'

        new File('/path/file_name').getBaseName() == 'file_name'
        Paths.get('/path/file_name').getBaseName() == 'file_name'

        new File('/path/').getBaseName() == 'path'
        Paths.get('/path/').getBaseName() == 'path'

        new File('/').getBaseName() == ''
        Paths.get('/').getBaseName() == ''

    }


    def testGetName() {

        expect:
        new File('/path/file.txt').getName() == 'file.txt'
        Paths.get('/path/file.txt').getName() == 'file.txt'

        new File('/path/file_name').getName() == 'file_name'
        Paths.get('/path/file_name').getName() == 'file_name'

        new File('/path/').getName() == 'path'
        Paths.get('/path/').getName() == 'path'

        new File('/').getName() == ''
        Paths.get('/').getName() == ''

    }

    def testPathMakeDir() {

        when:
        def path = Paths.get('xyz')
        then:
        path.mkdir()

        cleanup:
        path.deleteDir()


    }


    def 'test baseName' () {

        expect:
        new File('a/b/file.text').baseName == 'file'
        new File('.').baseName == ''

    }

    def 'test file extension' () {

        expect:
        new File('a/b/file.txt').extension == 'txt'
        new File('.').baseName == ''

    }

    def 'test copyTo' () {

        setup:
        def source = File.createTempFile('test','source')
        source.text = 'Hello'

        when:
        def copy = source.copyTo( File.createTempFile('test','target') )
        then:
        copy.text == 'Hello'
        source.exists()

        when:
        def folder = Files.createTempDirectory('testCopyTo').toFile()
        def copy2 = source.copyTo(folder)
        then:
        copy2.text == 'Hello'
        copy2.name == source.name
        source.exists()

        cleanup:
        source?.delete()
        copy?.delete()
        folder?.deleteDir()
    }

    def 'test copy source directory (file) to target directory' () {

        when:
        def sourceFolder =  Files.createTempDirectory('folderToCopy').toFile()
        def path1 = new File(sourceFolder, 'file1')
        def path2 = new File(sourceFolder, 'file2')
        def path3 = new File(sourceFolder, 'file3')
        path1.text = 'Hello1'
        path2.text = 'Hello2'
        path3.text = 'Hello3'
        def target = Files.createTempDirectory('targetFolder').resolve('xxx').toFile()
        // try to copy the source folder to the target
        sourceFolder.copyTo(target)

        then:
        new File(target, 'file1').text == 'Hello1'
        new File(target, 'file2').text == 'Hello2'
        new File(target, 'file3').text == 'Hello3'
        !new File(target, 'file4').exists()

        cleanup:
        sourceFolder?.deleteDir()
        target?.deleteDir()

    }


    def 'test copyTo path' () {

        setup:
        def source = Files.createTempFile('test','source')
        source.text = 'Hello'

        when:
        def copy = source.copyTo( Files.createTempFile('test','target') )
        then:
        copy.text == 'Hello'
        source.exists()

        when:
        def folder = Files.createTempDirectory('testCopyTo')
        def copy2 = source.copyTo(folder)
        then:
        copy2.text == 'Hello'
        copy2.name == source.name
        source.exists()

        cleanup:
        source?.delete()
        copy?.delete()
        folder?.deleteDir()
    }

    def 'test copy source directory (path) to another directory' () {

        when:
        def sourceFolder =  Files.createTempDirectory('folderToCopy')
        def path1 = sourceFolder.resolve('file1')
        def path2 = sourceFolder.resolve('file2')
        def path3 = sourceFolder.resolve('file3')
        def sub1 = Files.createTempDirectory(sourceFolder, 'sub1')
        def path4 = sub1.resolve('file4')
        path1.text = 'Hello1'
        path2.text = 'Hello2'
        path3.text = 'Hello3'
        path4.text = 'Hello4'
        def targetFolder = Files.createTempDirectory('targetFolder')
        // try to copy the source folder to the target
        sourceFolder.copyTo(targetFolder)
        then:
        targetFolder.resolve('file1').text == 'Hello1'
        targetFolder.resolve('file2').text == 'Hello2'
        targetFolder.resolve('file3').text == 'Hello3'
        targetFolder.resolve(sub1.getName()).resolve('file4').text == 'Hello4'

        cleanup:
        sourceFolder.deleteDir()
        targetFolder.deleteDir()

    }

    def 'test moveTo' () {

        setup:
        // source1
        def source1 = Files.createTempFile('test','source')
        source1.text = 'Hello 1'


        /*
         * test moving to a file with a different name
         */
        when:
        def file1 = source1.moveTo( Files.createTempFile('targetFile',null) )
        then:
        !source1.exists()
        file1.text == 'Hello 1'

        /*
         * test moving to a FOLDER
         */
        when:
        def folder = Files.createTempDirectory('testMoveTo')
        def source2 = Files.createTempFile('test','source')
        source2.text = 'Hello 2'
        def file2 = source2.moveTo(folder)
        then:
        file2.text == 'Hello 2'
        file2.name == source2.name
        !source2.exists()


        cleanup:
        file1?.delete()
        file2?.delete()
        folder?.deleteDir()

    }

    def 'test moveTo (directory) ' () {

        given:
        def sourceFolder =  Files.createTempDirectory('folderToCopy')
        def path1 = sourceFolder.resolve('file1')
        def path2 = sourceFolder.resolve('file2')
        def path3 = sourceFolder.resolve('file3')
        def sub1 = Files.createTempDirectory(sourceFolder, 'sub1')
        def path4 = sub1.resolve('file4')
        path1.text = 'Hello1'
        path2.text = 'Hello2'
        path3.text = 'Hello3'
        path4.text = 'Hello4'
        def targetFolder = Files.createTempDirectory('targetFolder')

        when:
        targetFolder = sourceFolder.moveTo(targetFolder)
        then:
        !sourceFolder.exists()
        targetFolder.getName() == sourceFolder.getName()
        targetFolder.resolve('file1').text == 'Hello1'
        targetFolder.resolve('file2').text == 'Hello2'
        targetFolder.resolve('file3').text == 'Hello3'
        targetFolder.resolve(sub1.getName()).resolve('file4').text == 'Hello4'

        cleanup:
        sourceFolder?.deleteDir()
        targetFolder.getParent()?.deleteDir()



    }



}
