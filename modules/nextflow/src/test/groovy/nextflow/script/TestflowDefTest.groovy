package nextflow.script

import nextflow.ast.NextflowDSL
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
import test.Dsl2Spec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestflowDefTest extends Dsl2Spec {

    static abstract class TestScript extends BaseScript {

        private injectSession() {
            try {
                def sess = binding.getSession()
                if( !sess )
                    return
                def f = this.class.superclass.superclass.getDeclaredField('session')
                f.setAccessible(true)
                f.set(this, sess)
            }
            catch (GroovyCastException e) {
                log.warn "Cant inject session - not a ScriptBinding context object"
            }
        }

        Object run() {
            injectSession()
            runScript()
            return this
        }

    }

    def 'should parse testflow' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            process foo {
                script:
                'echo hello'
            }
        
            testflow bar {
              given:
                x=1
              when:
                foo()
              then:
                true
            }
        '''

        when:
        def script = (TestScript)new GroovyShell(new ScriptBinding(), config).parse(SCRIPT).run()
        def meta = ScriptMeta.get(script)

        then:
        meta.definitions.size() == 2
        and:
        meta.getProcess('foo')
        and:
        meta.getTestflow('bar')
        meta.getTestflow('bar').@given !=null
        meta.getTestflow('bar').@when !=null
        meta.getTestflow('bar').@then !=null
        and:
        !meta.getTestflow('meh')
    }


}
