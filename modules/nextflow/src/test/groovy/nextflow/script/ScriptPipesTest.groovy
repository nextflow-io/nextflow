package nextflow.script

import java.nio.file.Files

import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner
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
            main: Channel.from('Hello') | map { it.reverse() } | (foo & bar)
            emit:
                foo.out
                bar.out
        }
        
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: Channel.from('Hola') | foo | map { it.reverse() } | bar
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: Channel.from('hello') | foo | bar 
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: Channel.from('hola') | foo | concat 
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: foo | map { it.reverse() }
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: foo | bar
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

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
            emit: Channel.from(1,2,3) | square | collect 
        }
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val.sort() == [1, 4, 9]
    }


    def 'should pipe branch output to concat operator' () {
        given:
        def SCRIPT ='''   
        Channel.from(10,20,30) | branch { foo: it <=10; bar: true } | concat 
        '''

        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result.val == 10
        result.val == 20
        result.val == 30
    }


    def 'should pipe branch output to processes' () {
        given:
        def SCRIPT ='''
        process foo {
          input: val x 
          input: val y
          output: val ret
          exec: ret=x*2+y
        }

        workflow {
           emit: Channel.from(10,20) | branch { foo: it <=10; bar: true } | foo 
        }
        '''

        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result.val == 40
    }

    def 'should compose custom funs' () {
        given:
        def SCRIPT = """
        process foo {
          output: val ret
          exec: ret=10
        }
        
        def bar(ch) {
          ch.map { it +1 }
        }

        workflow {
            emit: foo | bar | map{ it*2 } 
        }
        """

        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result.val == 22

    }

    def 'should compose custom funs/2' () {
        given:
        def SCRIPT = """
        process foo {
          output: 
            val x
            val y
          exec: 
            x=1; y=2
        }

        def bar(ch1, ch2) {
          ch1.combine(ch2)
        }

        workflow {
            emit: foo | bar | view
        }
        """

        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result.val == [1,2]

    }

    def 'should compose imported funs' () {
        given:
        def folder = Files.createTempDirectory('test')
        def MODULE = folder.resolve('module.nf')
        MODULE.text = '''
        process foo {
          output: val ret
          exec: ret=10
        }

        def bar(ch) {
          ch.map { it +1 }
        }

        '''
        def SCRIPT = """
        include{ foo; bar } from "$MODULE"

        workflow {
            emit: foo | bar | map{ it*3 }
        }
        """

        when:
        def result = new MockScriptRunner().setScript(SCRIPT).execute()
        then:
        result.val == 33


        cleanup:
        folder?.deleteDir()
    }

}
