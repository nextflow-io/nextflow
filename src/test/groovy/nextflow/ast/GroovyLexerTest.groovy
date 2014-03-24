/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.ast

import groovy.transform.InheritConstructors
import org.codehaus.groovy.antlr.GroovySourceToken
import spock.lang.Ignore
import spock.lang.Specification
import test.TestParser

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GroovyLexerTest extends Specification{

    @InheritConstructors
    static class Lexer extends CustomLexer {
        int replaceCount

        protected replace( GroovySourceToken token ) {
            replaceCount++
            super.replace(token)
        }
    }

    def disable() {

    }

    @Ignore
    def testLexer() {

        setup:
        def Lexer lexer
        def parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })

        def script = '''
            process hola {
               input: val x as y
               'echo hello'
            }
            '''

        when:
        parser.parse(script)
        then:
        lexer.replaceCount == 1

    }

    @Ignore
    def testLexerWithTwoProcess() {

        setup:
        def Lexer lexer
        def parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })

        def script = '''
            process hola {
               input: val x as y
               'echo hello'
            }

            process hola {
               input:
               val x
               val y as z
               output:
               file z as 'file.txt'

               'echo hello'
            }
            '''

        when:
        parser.parse(script)
        then:
        lexer.replaceCount == 3

    }


    @Ignore
    def testLexerWithinIf() {

        setup:
        def Lexer lexer
        def parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })

        def script = '''
            process hola {
               input: val x as y
               'echo hello'
            }

            if( 1 ) {

               process hola {
                   input:
                   val x
                   val y as z
                   output:
                   file z as 'file.txt'

               'echo hello'
               }

            }
            '''

        when:
        parser.parse(script)
        then:
        lexer.replaceCount == 3

    }

    @Ignore
    def testNoReplace() {

        setup:
        def Lexer lexer
        def parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })
        def script = '''
            process hola {
               'echo hello'
            }

            x as File
            '''

        when:
        parser.parse(script)
        then:
        lexer.replaceCount == 0


        /*
         * the 'as' in this input refers to the type
         */
        when:
        parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })
        script = '''
            process hola {
                input: file ( a as File )
               'echo hello'
            }

            x as File
            '''
        parser.parse(script)

        then:
        lexer.replaceCount == 0



        /*
         * the 'as' in this input refers to the type
         */
        when:
        parser = new TestParser(lexerFactory: { lexer=new Lexer(it) })
        script = '''
            process hola {
                input: file as 'x'
                exec:
                x =1
                x = file as File
                y = 2 as String
            }

            '''
        parser.parse(script)

        then:
        lexer.replaceCount == 1

    }




}
