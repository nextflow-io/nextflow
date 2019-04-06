/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CsvParserTest extends Specification {

    def 'should parse csv line' () {

        given:
        def parser = new CsvParser()
                    .setSeparator(SEP)
                    .setQuote(QUOTE)

        expect:
        parser.parse(LINE) == EXPECTED

        where:
        LINE                    | SEP   | QUOTE | EXPECTED
        'a,b,c'                 | null  | null  | ['a','b','c']
        'a,,c'                  | null  | null  | ['a',null,'c']
        'a,,'                   | null  | null  | ['a',null,null]
        '"1","2","3"'           | null  | null  | ['"1"','"2"','"3"']
        '"1,0","2.1","3,4"'     | null  | null  | ['"1','0"','"2.1"','"3','4"']
        'a,b,c'                 | null  | '"'   | ['a','b','c']
        '"1","2","3"'           | null  | '"'   | ['1','2','3']
        '"1,0","2.1","3,4"'     | null  | '"'   | ['1,0','2.1','3,4']
        'AxBxC'                 | 'x'   | null  | ['A','B','C']
    }

    def 'should return empty value' () {

        given:
        def parser = new CsvParser(empty: EMPTY)
                .setQuote(QUOTE)

        expect:
        parser.parse(LINE) == EXPECTED

        where:
        LINE        | EMPTY | QUOTE | EXPECTED
        'a,,c'      | null  | null  | ['a',null,'c']
        'a,,'       | null  | null  | ['a',null,null]
        ',,,c'      | null  | null  | [null,null,null,'c']

        'a,,'       | ''    | null  | ['a','','']
        'a,,'       | ''    | null  | ['a','','']
        ',,,c'      | ''    | null  | ['','','','c']

        'a,,'       | ''    | '"'   | ['a','','']
        'a,,'       | ''    | '"'   | ['a','','']
        ',,,c'      | ''    | '"'   | ['','','','c']
    }


    def 'should strip blanks'( ) {

        when:
        def splitter = new CsvParser().setQuote('"')
        then:
        splitter.stripBlanks('hello') == 'hello'
        splitter.stripBlanks('"hello"') == 'hello'
        splitter.stripBlanks('\'hello\'') == '\'hello\''
        splitter.stripBlanks('""') == ''
        splitter.stripBlanks('"') == '"'

        when:
        splitter = new CsvParser().setQuote('"').setStrip(true)
        then:
        splitter.stripBlanks('hello ') == 'hello'
        splitter.stripBlanks('"hello"') == 'hello'
        splitter.stripBlanks('" hello "  ') == ' hello '
        splitter.stripBlanks(' " " ') == ' '
        splitter.stripBlanks(' "" ') == ''

    }

}
