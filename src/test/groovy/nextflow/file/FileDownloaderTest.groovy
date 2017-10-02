package nextflow.file

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Session
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileDownloaderTest extends Specification {

    def 'should get the max threads value' () {

        when:
        new Session()
        then:
        FileDownloader.getMaxThreads() == Runtime.getRuntime().availableProcessors()

        when:
        new Session([download: [maxThreads: 99]])
        then:
        FileDownloader.getMaxThreads() == 99

    }

    def 'should get the max retries value' () {

        when:
        new Session()
        then:
        FileDownloader.getMaxRetries() == 3

        when:
        new Session([download: [maxRetries: 88]])
        then:
        FileDownloader.getMaxRetries() == 88

    }

    def 'should download foreign files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def session = new Session(); session.workDir = folder
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')
        def local = Paths.get('local.txt')
        def files = [foo: local, bar: foreign1, baz: foreign2]

        when:
        def d = new FileDownloader()
        def result = d.downloadForeignFiles(files)
        then:
        result.foo ==  Paths.get('local.txt')

        result.bar.name == 'hola.txt'
        result.bar.text == 'hola mundo!'
        result.bar.fileSystem == FileSystems.default

        result.baz.name == 'ciao.txt'
        result.baz.text == 'ciao mondo!'
        result.baz.fileSystem == FileSystems.default

        cleanup:
        folder?.deleteDir()
    }
}
