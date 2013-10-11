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

    def 'test moveTo' () {

        setup:
        def folder = Files.createTempDirectory('testMoveTo').toFile()
        // source1
        def source = File.createTempFile('test','source')
        source.text = 'Hello'
        // source2
        def source2 = File.createTempFile('test','source')
        source2.text = 'Hello'

        /*
         * test moving to a file with a different name
         */
        when:
        def file1 = source.moveTo( new File('target') )
        then:
        !source.exists()
        file1.text == 'Hello'
        file1.name == 'target'

        /*
         * test moving to a FOLDER
         */
        when:
        def file2 = source2.moveTo(folder)
        then:
        file2.text == 'Hello'
        file2.name == source2.name
        !source2.exists()


        cleanup:
        file1?.delete()
        file2?.delete()
        folder?.deleteDir()

    }


}
