package nextflow.script

import spock.lang.Timeout

import nextflow.Channel
import nextflow.exception.ScriptCompilationException
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ScriptDslTest extends Dsl2Spec {


    def 'should execute basic workflow' () {
        when:
        def result = dsl_eval '''
   
        workflow {
            emit: result 
            main:
            result = 'Hello world'
        }
        '''

        then:
        result.val == 'Hello world'
    }

    def 'should execute emit' () {
        when:
        def result = dsl_eval '''
   
        workflow {
            emit:
            result = 'Hello world'
        }
        '''

        then:
        result.val == 'Hello world'
    }

    def 'should emit expression' () {
        when:
        def result = dsl_eval '''
         
        def foo() { 'Hello world' } 
       
        workflow {
            emit:
            foo().toUpperCase()
        }
        '''

        then:
        result.val == 'HELLO WORLD'
    }

    def 'should emit process out' () {
        when:
        def result = dsl_eval '''
         
        process foo {
          output: val x 
          exec: x = 'Hello'
        }
       
        workflow {
            main: foo()
            emit: foo.out
        }
        '''

        then:
        result.val == 'Hello'
    }


    def 'should define processes and workflow' () {
        when:
        def result = dsl_eval '''
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
        
        workflow alpha {
            get: 
                data
            
            main:
                foo(data)
                bar(foo.out)
                
            emit: 
                x = bar.out 
            
        }
   
        workflow {
            main: alpha('Hello') 
            emit: x = alpha.out 
        }
        '''

        then:
        result.val == 'HELLO MUNDO'
    }


    def 'should access nextflow enabling property' () {
        when:
        def result = dsl_eval '''
        return nextflow.preview.dsl 
        '''

        then:
        result == 2
    }


    def 'should not allow function with reserved identifier' () {

        when:
        dsl_eval """ 
            def main() { println 'ciao' }
        """

        then:
        def err = thrown(ScriptCompilationException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow process with reserved identifier' () {

        when:
        dsl_eval """ 
            process main {
              /echo ciao/
            }
        """

        then:
        def err = thrown(ScriptCompilationException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow workflow with reserved identifier' () {

        when:
        dsl_eval """ 
            workflow main {
              /echo ciao/
            }
        """
        then:
        def err = thrown(ScriptCompilationException)
        err.message.contains('Identifier `main` is reserved for internal use')
    }

    def 'should not allow duplicate workflow keyword' () {
        when:
        dsl_eval(
                """ 
                workflow {
                  /echo ciao/
                }
                
                workflow {
                  /echo miao/
                }
            """
        )
        then:
        def err = thrown(ScriptCompilationException)
        err.message.contains('Duplicate entry workflow definition')
    }

    def 'should apply operator to process result' () {
        when:
        def result = dsl_eval(/
            process hello {
              output: val result
              exec:
                result = "Hello"
            }     
            
            workflow {
               main: hello()
               emit: hello.out.map { it.toUpperCase()  }
            }
        /)
        then:
        result.val == 'HELLO'
    }

    def 'should branch and view' () {

        when:
        def result = dsl_eval(/
            Channel
                .from(1,2,3,40,50)
                .branch { 
                    small: it < 10 
                    large: it > 10  
                }
                .set { result }
                
             ch1 = result.small.map { it }
             ch2 = result.large.map { it }  

             [ch1, ch2]
        /)
        then:
        result[0].val == 1
        result[0].val == 2
        result[0].val == 3
        and:
        result[1].val == 40
        result[1].val == 50
    }


    def 'should allow pipe process and operator' () {
        when:
        def result = dsl_eval('''
        process foo {
          output: val result
          exec: result = "hello"
        }     
 
        process bar {
          output: val result
          exec: result = "world"
        } 
        
        workflow {
           emit: (foo & bar) | concat      
        }
        ''')

        then:
        result.val == 'hello'
        result.val == 'world'
        result.val == Channel.STOP
    }

    def 'should allow process and operator composition' () {
        when:
        def result = dsl_eval('''
        process foo {
          output: val result
          exec: result = "hello"
        }     
 
        process bar {
          output: val result
          exec: result = "world"
        } 
        
        workflow {
           main: foo(); bar()
           emit: foo.out.concat(bar.out)      
        }
        ''')


        then:
        result.val == 'hello'
        result.val == 'world'
        result.val == Channel.STOP
    }

}
