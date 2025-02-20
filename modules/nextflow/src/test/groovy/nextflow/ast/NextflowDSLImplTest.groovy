package nextflow.ast

import nextflow.Session
import nextflow.script.BaseScript
import nextflow.script.ScriptMeta
import nextflow.script.ScriptParser
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import test.Dsl2Spec
import test.MockExecutorFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowDSLImplTest extends Dsl2Spec {

    def createCompilerConfig() {
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        return config
    }


    def 'should fetch method names' () {

        given:
        def config = createCompilerConfig()

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
        def config = createCompilerConfig()

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
        def config = createCompilerConfig()

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
        def config = createCompilerConfig()

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
        def session = new Session()
        session.executorFactory = new MockExecutorFactory()
        and:
        def parser = new ScriptParser(session)

        def SCRIPT = '''
            process alpha {
              /hello/
            }
        
            process beta {
              /world/
            }
            
            workflow {
              alpha(); beta()
            }
        '''

        when:
        parser.runScript(SCRIPT)
        then:
        ScriptMeta.get(parser.getScript()).getProcessNames() == ['alpha', 'beta'] as Set
    }

    def 'should fetch variable names' () {

        given:
        def config = createCompilerConfig()

        def SCRIPT = '''
            process alpha {
                input:
                val foo

                exec:
                if( foo == 'a' )
                    log.info "foo is 'a'"
                def bar = foo
            }

            workflow {
                alpha('a')
            }
            '''

        when:
        new GroovyShell(config).parse(SCRIPT)
        then:
        noExceptionThrown()
    }

}
