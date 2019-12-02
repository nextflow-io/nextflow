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

package nextflow.ast

import spock.lang.Specification

import nextflow.script.TokenBranchDef
import nextflow.script.TokenBranchChoice
import nextflow.script.TokenMultiMapDef
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class OpXformTest extends Specification {

    private TokenBranchDef eval_branch(String stmt) {

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( new ASTTransformationCustomizer(OpXform))

        def shell = new GroovyShell(config)
        def result = shell.evaluate("""
        class OpTest {
            def branch(Closure c) { c.call() }
        }
        
        new OpTest().branch ( $stmt )

        """)

        return (TokenBranchDef)result
    }

    private TokenMultiMapDef eval_multiMap(String stmt) {

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( new ASTTransformationCustomizer(OpXform))

        def shell = new GroovyShell(config)
        def result = shell.evaluate("""
        class OpTest {
            def multiMap(Closure c) { c.call() }
        }
        
        new OpTest().multiMap ( $stmt )

        """)

        return (TokenMultiMapDef)result
    }

    def 'should transform basic switch' () {
        when:
        def result = eval_branch('''{
            x: it < 1; return 'hello'
            y: it > 1; return 'world'
        }''')

        then:
        result instanceof TokenBranchDef
        result.closure instanceof Closure<TokenBranchChoice>
        result.branches == ['x','y']
        and:
        result.closure.call(0).value == 'hello'
        result.closure.call(0).choice == 'x'
        and:
        result.closure.call(2).value == 'world'
        result.closure.call(2).choice == 'y'
        and:
        result.closure.call(1) == null
    }

    def 'should transform switch statement' () {
        when:
        def result = eval_branch('''
            {
                x: it < 1; 'hello'
                y: it > 1; 'world'
                z: true; 'the-default-value'
            }
        ''')

        then:
        result.branches == ['x','y','z']
        result.closure(0) == new TokenBranchChoice('hello', 'x')
        result.closure(2) == new TokenBranchChoice('world', 'y')
        result.closure(1) == new TokenBranchChoice('the-default-value', 'z')
    }


    def 'should return implicit values' () {
        when:
        def result = eval_branch('''
            {
                x: it < 1
                y: it > 1
                z: true
            }
        ''')

        then:
        result.branches == ['x','y','z']
        result.closure(-10) == new TokenBranchChoice(-10, 'x')
        result.closure(20) == new TokenBranchChoice(20, 'y')
        result.closure(1) == new TokenBranchChoice(1, 'z')
    }

    def 'should return explicit param' () {
        when:
        def result = eval_branch('''
        { p ->
            x: p < 1
            y: p > 1
            z: true;
        }
        ''')

        then:
        result.branches == ['x','y','z']
        result.closure(-10) == new TokenBranchChoice(-10 , 'x')
        result.closure(20) == new TokenBranchChoice(20, 'y')
        result.closure(1) == new TokenBranchChoice(1, 'z')
    }


    def 'should allow arbitrary prolog code' () {
        when:
        def result = eval_branch('''
            { p ->
                def alpha=1
                def delta=2
                def omega=alpha+delta
                x: p < 1; 
                return omega;
            }
        ''')

        then:
        result.branches == ['x']
        result.closure.call(0).value == 3 // ie. omega = 1+2
        result.closure.call(0).choice == 'x'
    }


    def 'should parse fork block' () {
        when:
        def result = eval_multiMap('''
            { it -> 
                foo:
                it+1
                
                bar: 
                it*it
            }
        ''')

        then:
        result.names == ['foo', 'bar']
        result.closure.call(1) == [foo:2, bar:1]
        result.closure.call(2) == [foo:3, bar:4]
        result.closure.call(3) == [foo:4, bar:9]
    }


    def 'should parse fork multi-block' () {
        when:
        def result = eval_multiMap('''
            { it -> 
                alpha:
                beta:
                delta:
                it+1
                
                omega:
                it*2
            }
        ''')

        then:
        result.names == ['alpha', 'beta', 'delta', 'omega']
        result.closure.call(10) == [alpha:11, beta: 11, delta:11, omega:20]
    }


    def 'should parse fork long ending' () {
        when:
        def result = eval_multiMap('''
            { it -> 
                alpha:
                it+1
                omega:
                def x=1
                def y=2
                it+x+y
            }
        ''')

        then:
        result.names == ['alpha', 'omega']
        result.closure.call(10) == [alpha:11, omega:13]
    }

    def 'should parse empty fork' () {
        when:
        eval_multiMap('''
            { it -> it+1 }
        ''')

        then:
        def e = thrown(MultipleCompilationErrorsException)
        e.message.contains "The forking criteria should define at least two target channels"
    }
}
