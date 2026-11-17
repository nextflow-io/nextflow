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
package nextflow.script.formatter;

import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import static nextflow.script.formatter.FormatterTestUtils.assertFormat;
import static nextflow.script.formatter.FormatterTestUtils.assertFormatted;

/**
 * Specifies the output of {@link ScriptFormattingVisitor} for statements and
 * expressions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ScriptFormatterStatementTest {

    @Nested
    class Statements {

        @Test
        void formatAssert() {
            assertFormat(
                """
                assert 2+2==4:'The math broke!'
                """,
                """
                assert 2 + 2 == 4 : 'The math broke!'
                """
            );
        }

        @Test
        void formatVariableDeclarations() {
            assertFormat(
                """
                def x=42
                def(x,y)=tuple(1,2)
                def(
                x,
                y
                )=tuple(1,2)
                """,
                """
                def x = 42
                def (x, y) = tuple(1, 2)
                def (x, y) = tuple(1, 2)
                """
            );
        }

        @Test
        void convertLegacyVariableType() {
            assertFormat(
                """
                def Integer x=42
                """,
                """
                def x: Integer = 42
                """
            );
        }

        @Test
        void formatTypedVariableDeclarations() {
            assertFormat(
                """
                def x:Integer=42
                def(x:Integer,y:Integer)=tuple(1,2)
                """,
                """
                def x: Integer = 42
                def (x: Integer, y: Integer) = tuple(1, 2)
                """
            );
        }

        @Test
        void formatAssignments() {
            assertFormat(
                """
                v=42
                list[0]='first'
                map.key='value'
                (x,y)=tuple(1,2)
                (
                x,
                y
                )=tuple(1,2)
                """,
                """
                v = 42
                list[0] = 'first'
                map.key = 'value'
                (x, y) = tuple(1, 2)
                (x, y) = tuple(1, 2)
                """
            );
        }

        @Test
        void formatIfElse() {
            assertFormat(
                """
                if(x<0.5)println('You lost.')else if(x>0.5)println('You won!')else println('You tied?')
                """,
                """
                if (x < 0.5) {
                    println('You lost.')
                }
                else if (x > 0.5) {
                    println('You won!')
                }
                else {
                    println('You tied?')
                }
                """
            );
        }

        @Test
        void formatThrow() {
            assertFormatted("""
                throw new Exception('something failed!')
                """);
        }

        @Test
        void formatTypedCatch() {
            assertFormat(
                """
                try{println(file('foo.txt').text)}catch(IOException e){log.warn("Could not load foo.txt")}
                """,
                """
                try {
                    println(file('foo.txt').text)
                }
                catch (e: IOException) {
                    log.warn("Could not load foo.txt")
                }
                """
            );
        }

        @Test
        void formatUntypedCatch() {
            assertFormat(
                """
                try{println(file('foo.txt').text)}catch(e){log.warn("Could not load foo.txt")}
                """,
                """
                try {
                    println(file('foo.txt').text)
                }
                catch (e) {
                    log.warn("Could not load foo.txt")
                }
                """
            );
        }

        @Test
        void formatMethodChain() {
            assertFormat(
                """
                channel.of( 1, 2, 3 )
                    .multiMap{v->foo:bar:v}.set{result}
                """,
                """
                channel.of(1, 2, 3)
                    .multiMap { v -> foo: bar: v }
                    .set { result }
                """
            );
        }
    }

    @Nested
    class Expressions {

        @Test
        void formatNumericLiterals() {
            assertFormatted("""
                42
                -1
                0b1001
                031
                0xabcd
                3.14
                -0.1
                1.59e7
                1.59e-7
                """);
        }

        @Test
        void formatKeywordLiterals() {
            assertFormatted("""
                true
                false
                null
                """);
        }

        @Test
        void formatStringLiterals() {
            assertFormatted("""
                "I said 'hello'"

                'I said "hello" again!'

                '''
                Hello,
                How are you today?
                '''

                \"""
                We don't have to escape quotes anymore!
                Even "double" quotes!
                \"""

                /no escape!/
                """);
        }

        @Test
        void formatInterpolatedString() {
            assertFormatted("""
                "Hello, ${names.join(' and ')}!"
                """);
        }

        @Test
        void wrapBareInterpolatedProperties() {
            assertFormat(
                """
                "Hello, $name.first $name.last!"
                """,
                """
                "Hello, ${name.first} ${name.last}!"
                """
            );
        }

        @Test
        void preserveLineContinuationsInTripleQuotedString() {
            assertFormatted("""
                \"""
                blastp \\
                    -in ${input} \\
                    -out ${output} \\
                    -db ${blast_db} \\
                    -html
                \"""
                """);
        }

        @Test
        void preserveDollarBracesInSingleQuotedString() {
            assertFormatted("""
                'Hello, ${names.join(" and ")}!'
                """);
        }

        @Test
        void formatListLiterals() {
            assertFormat(
                """
                [1,2,3]
                []
                """,
                """
                [1, 2, 3]
                []
                """
            );
        }

        @Test
        void formatMapLiterals() {
            assertFormat(
                """
                [foo:1,bar:2,baz:3]
                [:]
                """,
                """
                [foo: 1, bar: 2, baz: 3]
                [:]
                """
            );
        }

        @Test
        void formatMapWithVariableKey() {
            assertFormat(
                """
                [(x):1]
                """,
                """
                [(x): 1]
                """
            );
        }

        @Test
        void formatMapWithPropertyKey() {
            assertFormat(
                """
                [(x.y):1]
                """,
                """
                [(x.y): 1]
                """
            );
        }

        @Test
        void formatClosure() {
            assertFormat(
                """
                {a,b->a+b}
                """,
                """
                { a, b -> a + b }
                """
            );
        }

        @Test
        void convertLegacyClosureParameterTypes() {
            assertFormat(
                """
                {Integer a,Integer b->a+b}
                """,
                """
                { a: Integer, b: Integer -> a + b }
                """
            );
        }

        @Test
        void formatTypedClosure() {
            assertFormat(
                """
                {a:Integer,b:Integer->a+b}
                """,
                """
                { a: Integer, b: Integer -> a + b }
                """
            );
        }

        @Test
        void formatMultiStatementClosure() {
            assertFormat(
                """
                {v->println'Hello!';v*v}
                """,
                """
                { v ->
                    println('Hello!')
                    v * v
                }
                """
            );
        }

        @Test
        void formatIndexAndPropertyAccesses() {
            assertFormatted("""
                myList[0]
                myFile.text
                myFiles*.text
                myFile?.text
                """);
        }

        @Test
        void addParenthesesToFunctionCalls() {
            assertFormat(
                """
                printf 'Hello %s!\\n', 'World'
                file 'hello.txt', checkIfExists: true
                """,
                """
                printf('Hello %s!\\n', 'World')
                file('hello.txt', checkIfExists: true)
                """
            );
        }

        @Test
        void formatMethodCallsWithClosureArgument() {
            assertFormat(
                """
                [1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }
                [1, 2, 3].each() { v -> println v }
                [1, 2, 3].each { v -> println v }
                """,
                """
                [1, 2, 3].inject('result:') { acc, v -> acc + ' ' + v }
                [1, 2, 3].each { v -> println(v) }
                [1, 2, 3].each { v -> println(v) }
                """
            );
        }

        @Test
        void formatConstructorCalls() {
            assertFormatted("""
                new java.util.Date()
                new Date()
                """);
        }

        @Test
        void formatUnaryBinaryAndTernaryExpressions() {
            assertFormat(
                """
                !(2+2==4)
                (1+2)*3
                x%2==0?'x is even!':'x is odd!'
                (false?'foo':true)?'bar':'baz'
                """,
                """
                !(2 + 2 == 4)
                (1 + 2) * 3
                x % 2 == 0 ? 'x is even!' : 'x is odd!'
                (false ? 'foo' : true) ? 'bar' : 'baz'
                """
            );
        }
    }

}
