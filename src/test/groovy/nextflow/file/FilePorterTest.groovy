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
class FilePorterTest extends Specification {

    def 'should get the max threads value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxThreads() == Runtime.getRuntime().availableProcessors()

        when:
        new Session([filePorter: [maxThreads: 99]])
        then:
        FilePorter.getMaxThreads() == 99

    }

    def 'should get the max retries value' () {

        when:
        new Session()
        then:
        FilePorter.getMaxRetries() == 3

        when:
        new Session([filePorter: [maxRetries: 88]])
        then:
        FilePorter.getMaxRetries() == 88

    }

    def 'should copy foreign files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def session = new Session(); session.workDir = folder
        def foreign1 = TestHelper.createInMemTempFile('hola.txt', 'hola mundo!')
        def foreign2 = TestHelper.createInMemTempFile('ciao.txt', 'ciao mondo!')
        def local = Paths.get('local.txt')
        def files = [foo: local, bar: foreign1, baz: foreign2]

        when:
        def d = new FilePorter()
        def result = d.stageForeignFiles(files)
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
