package nextflow.script

import spock.lang.Timeout
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ScriptPipesTest extends Dsl2Spec {

    def 'should pipe processes parallel' () {
        given:
        def SCRIPT =  '''
        process foo {
          input: val data 
          output: val result
          exec:
            result = "$data mundo"
        }     
        
        process bar {
            input: val data 
            output: val result
            exec: 
              result = data.toUpperCase()
        } 

        workflow {
            channel.of('Hello') | map { it.reverse() } | (foo & bar)
        }
        
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result[0].val == 'olleH mundo'
        result[1].val == 'OLLEH'
    }


    def 'should pipe processes sequential' () {
        given:
        def SCRIPT =  '''
        process foo {
          input: val data 
          output: val result
          exec:
            result = "$data world"
        }     
        
        process bar {
            input: val data 
            output: val result
            exec: 
              result = data.toUpperCase()
        } 

        workflow {
            channel.of('Hola') | foo | map { it.reverse() } | bar
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'DLROW ALOH'
    }

    def 'should pipe multi-outputs with multi-inputs' () {
        given:
        def SCRIPT =  '''
        process foo {
          input: val data 
          output: 
            val X
            val Y
          exec:
            X = data.reverse()
            Y = data.toUpperCase()
        }     
        
        process bar {
          input: 
            val X
            val Y 
          output: 
            val Z
          exec: 
            Z = "$X + $Y"  
        }
        
        // execute `foo` process and pipe
        // the multiple output channels 
        // to the `bar` process receiving multiple inputs
        workflow {
            channel.of('hello') | foo | bar 
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'olleh + HELLO'
    }

    def 'should pipe process with operator' () {
        given:
        def SCRIPT =  '''
        process foo {
          input: val data 
          output: 
            val X
            val Y
            val Z
          exec:
            X = data.reverse()
            Y = data.toUpperCase()
            Z = data
        }     
        
        // execute `foo` process and 
        // pipe the multiple output channels 
        // to the `concat` operator
        workflow {
            channel.of('hola') | foo | concat 
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'aloh'
        result.val == 'HOLA'
        result.val == 'hola'
    }


    def 'should allow process as first operand in a pipe chain' () {
        given:
        def SCRIPT =  '''
        process foo {
          output: 
            val X
          exec:
            X = "hola"
        }     
        
        workflow {
            foo | map { it.reverse() }
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'aloh'
    }


    def 'should allow two process piping' () {
        given:
        def SCRIPT =  '''
        process foo {
          output: val X
          exec:
            X = "hola"
        }     
        
        process bar {
          input: val X
          output: val Z
          exec:
            Z = X.toUpperCase()  
        }
        
        workflow {
            foo | bar
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'HOLA'
    }


    def 'should pipe with each input' () {
        given:
        def SCRIPT =  '''
        process square {
          input: each X
          output: val Z
          exec:
            Z = X*X
        }     
        
        workflow {
            channel.of(1,2,3) | square | collect 
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val.sort() == [1, 4, 9]
    }


    def 'should pipe branch output to concat operator' () {
        given:
        def SCRIPT ='''   
        channel.of(10,20,30) | branch { foo: it <=10; bar: true } | concat 
        '''

        when:
        def result = runScript(SCRIPT)
        then:
        result.val == 10
        result.val == 20
        result.val == 30
    }


    def 'should pipe branch output to processes' () {
        given:
        def SCRIPT ='''
        process foo {
          input: val x ; val y
          output: val ret
          exec: ret=x*2+y
        }

        workflow {
           channel.of(10,20) | branch { foo: it <=10; bar: true } | foo 
        }
        '''

        when:
        def result = runScript(SCRIPT)
        then:
        result.val == 40
    }

}
