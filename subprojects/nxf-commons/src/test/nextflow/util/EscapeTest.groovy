package nextflow.util

import java.nio.file.Paths

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EscapeTest extends Specification {

    def 'should escape quotes in file names' () {
        expect:
        Escape.path(Paths.get('hello.txt')) == "hello.txt"
        Escape.path(Paths.get("hello'3.txt")) == "hello\\'3.txt"
        Escape.path(Paths.get("hello'3.txt")).size() == "hello'3.txt".size()+1
        Escape.path(Paths.get("/some'5/data'3/with/quote's/file's.txt")) == "/some\\'5/data\\'3/with/quote\\'s/file\\'s.txt"
    }

    def 'should escape quote in file names as string' () {
        given:
        String world = 'world'

        expect:
        Escape.path('hello.txt') == "hello.txt"
        Escape.path("hello'3.txt") == "hello\\'3.txt"
        Escape.path("hello'3.txt").size() == "hello'3.txt".size()+1
        Escape.path("/some'5/data'3/with/quote's/file's.txt") == "/some\\'5/data\\'3/with/quote\\'s/file\\'s.txt"
        Escape.path("Hello '$world'") == "Hello\\ \\'world\\'"

    }

}
