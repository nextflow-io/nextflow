/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

}
