package nextflow.script

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.ast.NextflowDSL
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
import spock.lang.Specification
import spock.lang.Timeout
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Timeout(5)
class WorkflowDefTest extends Specification {

    def setupSpec() { NextflowMeta.instance.enableDsl2() }
    def cleanupSpec() { NextflowMeta.instance.disableDsl2() }

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

    def 'should parse workflow' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              print 'Hello world'
            }
        
            workflow bravo(foo, bar) {
              print foo
              print bar
              return foo+bar
            }
            
            workflow delta(foo) {
                println foo+bar
            }

            workflow empty { }
        '''

        when:
        def script = (TestScript)new GroovyShell(config).parse(SCRIPT).run()
        def meta = ScriptMeta.get(script)
        then:
        meta.definitions.size() == 4
        meta.getWorkflow('alpha') .declaredInputs == []
        meta.getWorkflow('alpha') .declaredVariables == []
        meta.getWorkflow('alpha') .source.stripIndent() == "print 'Hello world'\n"

        meta.getWorkflow('bravo') .declaredInputs == ['foo', 'bar']
        meta.getWorkflow('bravo') .declaredVariables == []
        meta.getWorkflow('bravo') .source.stripIndent() == "print foo\nprint bar\nreturn foo+bar\n"

        meta.getWorkflow('delta') .declaredInputs == ['foo']
        meta.getWorkflow('delta') .declaredVariables == ['bar']

        meta.getWorkflow('empty') .source == ''
        meta.getWorkflow('empty') .declaredInputs == []
        meta.getWorkflow('empty') .declaredVariables == []
    }

    def 'should define anonymous workflow' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow {
              print 1
              print 2
            }
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (TestScript)new GroovyShell(binding, config).parse(SCRIPT).run()
        def meta = ScriptMeta.get(script)
        then:
        meta.getWorkflow(null).getSource().stripIndent() == 'print 1\nprint 2\n'

    }

    def 'should run workflow block' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha(x) {
              return "$x world"
            }
       
        '''

        when:
        def script = (TestScript)new GroovyShell(config).parse(SCRIPT).run()
        def workflow = ScriptMeta.get(script).getWorkflow('alpha')
        then:
        workflow.declaredInputs == ['x']

    }

    def 'should not fail' () {
        given:
        // this test checks that the closure used to define the workflow
        // does NOT define an implicit `it` parameter that would clash
        // with the `it` used by the inner closure

        def SCRIPT = """       
        
        workflow {
            Channel.empty().map { id -> id +1 }  
            Channel.empty().map { it -> def id = it+1 }  
        }
        """

        when:
        new MockScriptRunner().setScript(SCRIPT).execute()

        then:
        noExceptionThrown()
    }

    def 'should validate collect output'() {

        given:
        def workflow = new WorkflowDef(Mock(BaseScript), Mock(TaskBody))
        def result

        when:
        result = workflow.collectOutputs(null)
        then:
        result instanceof DataflowVariable
        result.val == null

        when:
        def array = new ChannelArrayList()
        result = workflow.collectOutputs(array)
        then:
        result == array

        when:
        def dataVar = new DataflowVariable()
        result = workflow.collectOutputs(dataVar)
        then:
        result == dataVar

        when:
        def dataQueue = new DataflowQueue()
        result = workflow.collectOutputs(dataQueue)
        then:
        result == dataQueue

        when:
        def dataBroad = new DataflowBroadcast()
        result = workflow.collectOutputs(dataBroad)
        then:
        result == dataBroad

        when:
        result = workflow.collectOutputs(['a', 'b'])
        then:
        result.size() == 2 
        result[0] instanceof DataflowVariable
        result[0].val == 'a'
        result[1] instanceof DataflowVariable
        result[1].val == 'b'

        when:
        result = workflow.collectOutputs(['a', dataQueue])
        then:
        result.size() == 2
        result[0] instanceof DataflowQueue
        result[0].val == 'a'
        result[1] instanceof DataflowQueue

        when:
        def var1 = new DataflowQueue(); var1 << 'a'
        def var2 = new DataflowQueue(); var2 << 'b'
        def var3 = new DataflowQueue(); var3 << 'c'
        result = workflow.collectOutputs([var1, new ChannelArrayList([var2, var3])])
        then:
        result.size() == 3
        result[0].val == 'a'
        result[1].val == 'b'
        result[2].val == 'c'
    }

    def 'should clone with a new name' () {
        given:
        def work = new WorkflowDef(Mock(BaseScript), Mock(TaskBody), 'woo')

        when:
        def copy = work.withName('bar')
        then:
        copy.getName() == 'bar'
    }

}
