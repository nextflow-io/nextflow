package nextflow.ast


import nextflow.script.BaseScript
import nextflow.script.ScriptMeta
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import test.BaseSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowDSLImplTest extends BaseSpec {

    def 'should fetch method names' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
            def foo() { 
                return 0 
            }
            
            def bar() { 
                return 1 
            }
            
            private baz() { 
                return 2 
            }
            
            process alpha {
              /hello/
            }

        '''
        when:
        new GroovyShell(config).parse(SCRIPT)
        then:
        noExceptionThrown()
    }

    def 'should not allow duplicate processes' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            process alpha {
              /hello/
            }
        
            process alpha {
              /world/
            }

        '''

        when:
        new GroovyShell(config).parse(SCRIPT)
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains 'Identifier `alpha` is already used by another definition'
    }


    def 'should not allow duplicate workflow' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              /hello/
            }
        
            workflow alpha {
              /world/
            }

        '''

        when:
        new GroovyShell(config).parse(SCRIPT)
        then:
        def err = thrown(MultipleCompilationErrorsException)
        err.message.contains 'Identifier `alpha` is already used by another definition'
    }

    def 'should throw illegal name exception' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        when:
        def SCRIPT = '''
            workflow 'foo:bar' {
              /hello/
            }

        '''
        new GroovyShell(config).parse(SCRIPT)
        then:
        def e = thrown(MultipleCompilationErrorsException)
        e.message.contains 'Process and workflow names cannot contain colon character'


        when:
        SCRIPT = '''
            process 'foo:bar' {
              /hello/
            }

        '''
        new GroovyShell(config).parse(SCRIPT)
        then:
        e = thrown(MultipleCompilationErrorsException)
        e.message.contains 'Process and workflow names cannot contain colon character'
    }

    def 'should set process name in the script meta' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            process alpha {
              /hello/
            }
        
            process beta {
              /world/
            }

        '''

        when:
        def script = new GroovyShell(config).parse(SCRIPT)
        then:
        ScriptMeta.get(script).getDsl1ProcessNames() == ['alpha', 'beta']
        ScriptMeta.get(script).getProcessNames() == ['alpha', 'beta'] as Set
    }



}
