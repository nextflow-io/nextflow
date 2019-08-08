package nextflow.script

import spock.lang.Specification

import nextflow.NextflowMeta
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptPipesTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }
    
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

        
        Channel.from('Hello') | map { it.reverse() } | (foo & bar)
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

        
        Channel.from('Hola') | foo | map { it.reverse() } | bar
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
        Channel.from('hello') | foo | bar 
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
        Channel.from('hola') | foo | concat 
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
        
        foo | map { it.reverse() }
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
        
        foo | bar
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
        
        Channel.from(1,2,3) | square | collect 
        '''

        when:
        def runner = new MockScriptRunner()
        def result = runner.setScript(SCRIPT).execute()

        then:
        result.val.sort() == [1, 4, 9]

    }

    def 'should set a channel in the global context' () {
        given:
        def SCRIPT =  '''
        Channel.from(1,2,3) | set { foo } 
        Channel.value('hola') | set { bar } 
        '''

        when:
        def runner = new MockScriptRunner()
        runner.setScript(SCRIPT).execute()

        then:
        runner.session.binding.foo.val == 1
        runner.session.binding.foo.val == 2
        runner.session.binding.foo.val == 3
        runner.session.binding.bar.value == 'hola'

    }

}
