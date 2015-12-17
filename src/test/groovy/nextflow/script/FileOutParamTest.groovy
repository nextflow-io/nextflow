package nextflow.script

import java.nio.file.Paths

import nextflow.exception.IllegalFileException
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileOutParamTest extends Specification {

    def 'should return a name relative to the workDir' () {

        given:
        def workDir = Paths.get('/a/b/c')

        expect:
        FileOutParam.relativizeName( Paths.get('hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativizeName( Paths.get('sub/dir/hello.txt'), workDir ) == 'sub/dir/hello.txt'
        FileOutParam.relativizeName( Paths.get('/a/b/c/hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativizeName( Paths.get('/a/b/c/some/path/hello.txt'), workDir ) == 'some/path/hello.txt'

        when:
        FileOutParam.relativizeName( Paths.get('/c/b/a/hello.txt'), workDir )
        then:
        thrown(IllegalFileException)
    }

    def 'should return a relative path' () {
        expect:
        FileOutParam.clean('hola.txt') == 'hola.txt'
        FileOutParam.clean('/hola.txt') == 'hola.txt'
        FileOutParam.clean('///hola.txt') == 'hola.txt'
        FileOutParam.clean('/hola/world.txt') == 'hola/world.txt'
    }

}
