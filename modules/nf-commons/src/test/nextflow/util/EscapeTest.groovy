/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.util

import java.nio.file.Paths

import spock.lang.Specification
import spock.lang.Unroll

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
        Escape.path(Paths.get("hello(3).txt")) == "hello\\(3\\).txt"
        Escape.path(Paths.get("hello\\3.txt")) == "hello\\\\3.txt"
        Escape.path(Paths.get("/some'5/data'3/with/quote's/file's.txt")) == "/some\\'5/data\\'3/with/quote\\'s/file\\'s.txt"
    }

    def 'should escape quote in file names as string' () {
        given:
        String world = 'world'

        expect:
        Escape.path('hello.txt') == "hello.txt"
        Escape.path("hello'3.txt") == "hello\\'3.txt"
        Escape.path("hello'3.txt").size() == "hello'3.txt".size()+1
        Escape.path("hello(3).txt") == "hello\\(3\\).txt"
        Escape.path("hello!3.txt") == "hello\\!3.txt"
        Escape.path('hello[!x].txt') == 'hello[!x].txt' // <-- this `!` should not be escaped because it's a glob negation http://man7.org/linux/man-pages/man7/glob.7.html
        Escape.path('hello[x!].txt') == 'hello[x\\!].txt'
        Escape.path('hello[!x.txt') == 'hello[\\!x.txt'
        Escape.path('hello![x,*.txt]') == 'hello\\![x,*.txt]'
        Escape.path("hello&3.txt") == "hello\\&3.txt"
        Escape.path("hello<3.txt") == "hello\\<3.txt"
        Escape.path("hello>3.txt") == "hello\\>3.txt"
        Escape.path("hello`3.txt") == "hello\\`3.txt"
        Escape.path('hello:3.txt') == "hello\\:3.txt"
        Escape.path('hello;3.txt') == "hello\\;3.txt"
        Escape.path('hello$3.txt') == 'hello\\$3.txt'
        Escape.path("/some'5/data'3/with/quote's/file's.txt") == "/some\\'5/data\\'3/with/quote\\'s/file\\'s.txt"
        Escape.path("Hello '$world'") == "Hello\\ \\'world\\'"

    }

    def 'should escape wildcards' () {

        expect:
        Escape.wildcards('file_*') == 'file_\\*'
        Escape.wildcards('file_??') == 'file_\\?\\?'
        Escape.wildcards('file_{a,b}') == 'file_\\{a,b\\}'
        Escape.wildcards('file_!a.txt') == 'file_\\!a.txt'
        Escape.wildcards('file_[!a].txt') == 'file_\\[\\!a\\].txt'
    }

    def 'should escape cli' () {
        expect:
        Escape.cli('nextflow','run','this') == 'nextflow run this'
        Escape.cli('nextflow','--foo','file.txt') == 'nextflow --foo file.txt'
        Escape.cli('nextflow','--foo','*.txt') == "nextflow --foo '*.txt'"
        Escape.cli('nextflow','--foo','*_{1,2}.fq') == "nextflow --foo '*_{1,2}.fq'"
        Escape.cli('nextflow','--foo','a b c') == "nextflow --foo 'a b c'"
        Escape.cli('nextflow','--foo','ab\'c') == "nextflow --foo ab\\'c"
        Escape.cli('nextflow','--foo','a[b-c]') == "nextflow --foo 'a[b-c]'"
        Escape.cli('nextflow','--foo','a[!b-c]') == "nextflow --foo 'a[!b-c]'"
    }

    @Unroll
    def 'should escape blanks' () {
        expect:
        Escape.blanks(STR) == EXPECT
        where:
        STR         | EXPECT
        'foo '      | 'foo '
        'foo\n'     | 'foo\\n'
        'foo\t'     | 'foo\\t'
        'foo\f'     | 'foo\\f'
        'foo\r'     | 'foo\\r'
    }

    def 'should escape special char' () {
        expect:
        Escape.variable(STR) == EXPECT
        where:
        STR         | EXPECT
        'foo'       | 'foo'
        'foo[x]bar' | 'foo[x]bar'
        'foo '      | 'foo '
        'foo:bar'   | 'foo:bar'
        'foo!bar'   | 'foo!bar'
        'foo[!x]bar'| 'foo[!x]bar'
        'foo|bar'   | 'foo|bar'
        'foo`bar'   | 'foo`bar'
        'foo&bar'   | 'foo&bar'
        'foo(x)bar' | 'foo(x)bar'
        'foo<x>bar' | 'foo<x>bar'
        and:
        '$foo'      | '\\$foo'
    }

    def 'should escape newlines' () {
        expect:
        Escape.path('a\nb') == 'a$\'\\n\'b'
        Escape.wildcards('a\nb') == 'a$\'\\n\'b'
    }

    def 'should escape double-quoted string' () {
        expect:
        Escape.doubleQuoted('a\\b"c$d`e f') == 'a\\\\b\\"c\\$d\\`e f'
    }

    @Unroll
    def 'should preserve a file name with special characters in bash' () {
        given:
        def name = 'a\nb c;d `e` $(f) "g" \'h ${i=$(j)}${i@P} $k'

        when:
        def proc = ['bash', '-c', CMD.replace('NAME', ARG(name))].execute()

        then:
        proc.text == name
        proc.waitFor() == 0

        where:
        CMD                     | ARG
        'printf %s NAME'        | { String it -> Escape.path(it) }
        'eval "printf %s NAME"' | { String it -> Escape.doubleQuoted(Escape.path(it)) }
    }

    @Unroll
    def 'should match output pattern #PATTERN in eval as the output file collector does' () {
        given:
        def folder = java.nio.file.Files.createTempDirectory('test')
        (OTHERS + FILE).each { folder.resolve(it).text = 'x' }
        def splitter = nextflow.file.FilePatternSplitter.glob().parse(PATTERN)
        def collected = { String name ->
            splitter.isPattern()
                ? java.nio.file.FileSystems.default.getPathMatcher('glob:' + PATTERN).matches(Paths.get(name))
                : splitter.strip(PATTERN) == name
        }

        when:
        def script = 'eval "ls -1d ' + Escape.doubleQuoted(Escape.pattern(PATTERN)) + '"'
        def proc = new ProcessBuilder('bash', '-c', script).directory(folder.toFile()).start()

        then:
        proc.inputStream.text == FILE + '\n'
        proc.waitFor() == 0
        and:
        collected(FILE)
        !OTHERS.any(collected)

        cleanup:
        folder?.deleteDir()

        where:
        PATTERN             | FILE                  | OTHERS
        'g\\[1\\].txt'      | 'g[1].txt'            | ['g1.txt']
        '\\{a,b\\}.txt'     | '{a,b}.txt'           | ['a.txt', 'b.txt']
        'f\\*.txt'          | 'f*.txt'              | ['f1.txt']
        'f\\?.txt'          | 'f?.txt'              | ['f1.txt']
        'x\\\\y.txt'        | 'x\\\\y.txt'          | ['x\\y.txt']
        'x\\:y\\!z.txt'     | 'x\\:y\\!z.txt'       | ['x:y!z.txt']
        'd$x.txt'           | 'd$x.txt'             | ['d.txt']
        'd$*'               | 'd$x.txt'             | ['dx.txt']
        '*.csv'             | 'data.csv'            | ['data.txt']
        '*\\.txt'           | 'a.txt'               | ['a.csv']
        'a\\b*'             | 'ab.txt'              | ['a\\b.txt']
        'x\\\\y*'           | 'x\\y.txt'            | ['xy.txt']
        'a b;c `d` "e".txt' | 'a b;c `d` "e".txt'   | ['a']
        'n\nm.txt'          | 'n\nm.txt'            | ['n', 'm.txt']
    }
}
