/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import spock.lang.Specification

import java.nio.file.Paths

import nextflow.processor.TaskPath
import nextflow.script.BaseScript
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowXformImplTest extends Specification {


    def 'should fix equals for paths' () {

        given:
        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))

        final p = Paths.get('/foo/bar/data.txt')
        final t = new TaskPath(p)
        final u = new TaskPath(p)
        final q = new TaskPath(p, 'hola.txt')
        final s = new TaskPath(p, 'ciao.txt')
        final z = new TaskPath(p, 'ciao.txt')

        expect:
        // although `equals` is true the `==` operator returns false
        // because `Path` implements `Comparable` interface, therefore
        // groovy uses `compareTo` instead of `equals`
        // See
        //   DefaultTypeTransformation#compareEqual(Object,Object)
        //   https://stackoverflow.com/questions/28355773/in-groovy-why-does-the-behaviour-of-change-for-interfaces-extending-compar#comment45123447_28387391
        t != p
        p != t

        // AST transformation fix the equality issue
        when:
        def shell = new GroovyShell(new Binding(p:p, t:t, q:q, s:s, u:u, z:z), config)
        shell.evaluate('''
            assert p == t && t == p
            assert t == u && u == t && t == t 
            assert s == z && z == s
            assert p != q && q != p
            assert t != q && q != t   
            assert q != s
            assert p != q 
            assert p != s  
        ''')
        then:
        noExceptionThrown()

    }



    def 'should fix declaration expression' () {

        given:
        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
        def shell = new GroovyShell(config)

        // AST transformation fix the equality issue
        when:
        shell.evaluate('''
            import java.nio.file.Path
            import java.nio.file.Paths
            import nextflow.processor.TaskPath

            def p = Paths.get('/foo/x')
            def t = new TaskPath(p)
            final x = (t == new TaskPath(p))
            assert x

        ''')
        then:
        noExceptionThrown()

        when:
        shell.evaluate('''
            import java.nio.file.Path
            import java.nio.file.Paths
            import nextflow.processor.TaskPath

            def p = Paths.get('/foo/x')

            final a = new TaskPath(p)
            final b = new TaskPath(p)
            final c = new TaskPath(p)
            final x = ( a == b && b == c )  
            assert x
  
        ''')
        then:
        noExceptionThrown()

    }

    def 'should compare mem units' () {

        given:
        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
        def shell = new GroovyShell(config)

        // AST transformation fix the equality issue
        when:
        shell.evaluate('''

           // memory as LEFT operand   
            assert 1.KB < 2000
            assert 1.KB < 2.KB
            assert 1.KB < '2 KB'

            assert 1.KB <= 1024
            assert 1.KB <= 1.KB
            assert 1.KB <= '1 KB'

            assert 1.KB >= 1024
            assert 1.KB >= 1.KB
            assert 1.KB >= '1 KB'

            assert 1.KB > 100
            assert 1.KB > 100.B
            assert 1.KB > '100 B'

            assert 1.KB == 1024
            assert 1.KB == 1.KB
            assert 1.KB == '1 KB'

            assert 1.KB != 2000
            assert 1.KB != 2.KB
            assert 1.KB != '2 KB'

            // memory as RIGHT operand   
            assert 100 < 1.GB 
            assert 100.B < 1.GB
            assert '100B' < 1.GB

            assert 100 <= 1.GB 
            assert 100.B <= 1.GB
            assert '100B' <= 1.GB

            assert 100 <= 100.GB 
            assert 100.B <= 100.GB
            assert '100B' <= 100.GB

            assert 100 >= 100.B 
            assert 100.B >= 100.B
            assert '100B' >= 100.B

            assert 200 > 100.B 
            assert 200.B > 100.B
            assert '200 B' > 100.B 
            
            assert 100 == 100.B
            assert 100.B == 100.B 
            assert '100 KB' == 100.KB

            assert 101 != 100.B
            assert 101.B != 100.B 
            assert '101 KB' != 100.KB
            
        ''')
        then:
        noExceptionThrown()

    }

    def 'should compare duration units' () {

        given:
        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
        def shell = new GroovyShell(config)

        // AST transformation fix the equality issue
        when:
        shell.evaluate('''

           // duration as LEFT operand   
            assert 1.sec < 10_000
            assert 1.sec < 10.sec
            assert 1.sec < '10 sec'

            assert 1.sec <= 1_000
            assert 1.sec <= 1.sec
            assert 1.sec <= '1 sec'

            assert 1.sec >= 1_000
            assert 1.sec >= 1.sec
            assert 1.sec >= '1 sec'

            assert 1.min >= 1_000
            assert 1.min >= 1.sec
            assert 1.min >= '1 sec'
  
            assert 1.sec == 1_000
            assert 1.sec == 1.sec
            assert 1.sec == '1 sec'
 
            assert 1.sec != 2_000
            assert 1.sec != 2.sec
            assert 1.sec != '2 sec'
            
            // duration as RIGHT operand   
            assert 1000 < 10.sec 
            assert 1.sec < 10.sec 
            assert '1 sec' < 10.sec 
            
            assert 10000 <= 10.sec 
            assert 10.sec <= 10.sec 
            assert '10 sec' <= 10.sec 

            assert 10000 >= 10.sec 
            assert 10.sec >= 10.sec 
            assert '10 sec' >= 10.sec 

            assert 10000 == 10.sec 
            assert 10.sec == 10.sec 
            assert '10 sec' == 10.sec 

            assert 11000 != 10.sec 
            assert 11.sec != 10.sec 
            assert '11 sec' != 10.sec                                                                      
        ''')
        then:
        noExceptionThrown()

    }

}
