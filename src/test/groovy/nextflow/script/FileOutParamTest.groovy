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
        FileOutParam.relativize( Paths.get('x'), workDir ) == 'x'
        FileOutParam.relativize( Paths.get('hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativize( Paths.get('sub/dir/hello.txt'), workDir ) == 'sub/dir/hello.txt'
        FileOutParam.relativize( Paths.get('/a/b/c/x'), workDir ) == 'x'
        FileOutParam.relativize( Paths.get('/a/b/c/hello.txt'), workDir ) == 'hello.txt'
        FileOutParam.relativize( Paths.get('/a/b/c/some/path/hello.txt'), workDir ) == 'some/path/hello.txt'

        when:
        FileOutParam.relativize( Paths.get('/c/b/a/hello.txt'), workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( Paths.get('/a/b/c'), workDir )
        then:
        thrown(IllegalFileException)
    }


    def 'should return a name relative to the workDir (with string)' () {

        given:
        def workDir = Paths.get('/a/b/c/')

        expect:
        FileOutParam.relativize( 'x', workDir ) == 'x'
        FileOutParam.relativize( 'hello.txt', workDir ) == 'hello.txt'
        FileOutParam.relativize( 'sub/dir/hello.txt', workDir ) == 'sub/dir/hello.txt'
        FileOutParam.relativize( '/a/b/c/x', workDir ) == 'x'
        FileOutParam.relativize( '/a/b/c/hello.txt', workDir ) == 'hello.txt'
        FileOutParam.relativize( '/a/b/c/some/path/hello.txt', workDir ) == 'some/path/hello.txt'

        when:
        FileOutParam.relativize( '/c/b/a/hello.txt', workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( '/a/b/c', workDir )
        then:
        thrown(IllegalFileException)

        when:
        FileOutParam.relativize( '/a/b/c/', workDir )
        then:
        thrown(IllegalFileException)

    }


}
