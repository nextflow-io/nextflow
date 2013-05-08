package nextflow.extension

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileExtensionsTest extends Specification {

    def 'test empty' () {

        setup:
        def fileEmpty = File.createTempFile('test','test')
        fileEmpty.deleteOnExit()

        def folderEmpty = File.createTempDir()

        def folderNotEmpty = File.createTempDir()
        def fileInFolder = new File(folderNotEmpty, 'filename')
        fileInFolder.createNewFile()

        def fileNotEmpty = File.createTempFile('test','test')
        fileNotEmpty.text = 'Hola'
        fileNotEmpty.deleteOnExit()

        expect:
        new File('non existing').isEmpty()
        fileEmpty.isEmpty()
        !fileNotEmpty.isEmpty()
        folderEmpty.isEmpty()
        !folderNotEmpty.isEmpty()

        cleanup:
        folderNotEmpty.deleteDir()
        folderEmpty.deleteDir()

    }


    def 'test notEmpty' () {

        setup:
        def fileEmpty = File.createTempFile('test','test')
        fileEmpty.deleteOnExit()

        def folderEmpty = File.createTempDir()

        def folderNotEmpty = File.createTempDir()
        def fileInFolder = new File(folderNotEmpty, 'filename')
        fileInFolder.createNewFile()

        def fileNotEmpty = File.createTempFile('test','test')
        fileNotEmpty.text = 'Hola'
        fileNotEmpty.deleteOnExit()

        expect:
        ! new File('non existing').isNotEmpty()
        !fileEmpty.isNotEmpty()
        fileNotEmpty.isNotEmpty()

        !folderEmpty.isNotEmpty()
        folderNotEmpty.isNotEmpty()

        cleanup:
        folderNotEmpty.deleteDir()
        folderEmpty.deleteDir()

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
