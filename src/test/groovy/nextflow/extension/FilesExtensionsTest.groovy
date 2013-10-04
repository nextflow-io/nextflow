package nextflow.extension

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
        def copy = File.createTempFile('test','target')

        when:
        source.copyTo( copy )

        then:
        copy.text == 'Hello'

        cleanup:
        source.delete()
        copy.delete()

    }

    def 'test moveTo' () {

        setup:
        def source = File.createTempFile('test','source')
        source.text = 'Hello'

        when:
        def copy = new File('target')
        source.moveTo( copy )

        then:
        !source.exists()
        copy.text == 'Hello'

        cleanup:
        copy.delete()

    }


}
