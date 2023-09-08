package nextflow.ast

import nextflow.Session
import nextflow.file.FileHelper
import nextflow.script.BaseScript
import nextflow.script.ScriptMeta
import nextflow.script.ScriptParser
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import spock.lang.Unroll
import test.Dsl2Spec

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowDSLImplTest extends Dsl2Spec {

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

    @Unroll
    def 'should throw illegal name exception when a variable has the same name as the script file' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        def scriptPath = FileHelper.getLocalTempPath().resolve('hello.nf')

        when:
        scriptPath.text = '''
            hello = 1
            println hello
            '''
        new GroovyShell(config).parse(scriptPath.toFile())
        then:
        def e = thrown(MultipleCompilationErrorsException)
        e.message.contains 'Cannot declare a variable identifier with the same name as the script file'

        when:
        scriptPath.text = '''
            def hello() {
                hello = 1
                println hello
            }
            
            hello()
            '''
        new GroovyShell(config).parse(scriptPath.toFile())
        then:
        e = thrown(MultipleCompilationErrorsException)
        e.message.contains 'Cannot declare a variable identifier with the same name as the script file'

        when:
        scriptPath.text = '''
            def hello = 1
            println hello
            '''
        new GroovyShell(config).parse(scriptPath.toFile())
        then:
        noExceptionThrown()

        when:
        scriptPath.text = '''
            def hello() {
                def hello = 1
                println hello
            }
            
            hello()
            '''
        new GroovyShell(config).parse(scriptPath.toFile())
        then:
        noExceptionThrown()

        cleanup:
        scriptPath.delete()
    }

    def 'should set process name in the script meta' () {
        given:
        def parser = new ScriptParser(new Session())

        def SCRIPT = '''
            process alpha {
              /hello/
            }
        
            process beta {
              /world/
            }
            
            workflow {}
            '''

        when:
        parser.runScript(SCRIPT)
        then:
        ScriptMeta.get(parser.getScript()).getProcessNames() == ['alpha', 'beta'] as Set
    }

}
