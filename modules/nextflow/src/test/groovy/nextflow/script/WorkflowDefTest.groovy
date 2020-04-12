package nextflow.script

import spock.lang.Timeout

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Session
import nextflow.ast.NextflowDSL
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
import test.Dsl2Spec
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Timeout(5)
class WorkflowDefTest extends Dsl2Spec {

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
        
            workflow bravo() {
              take: foo
              take: bar
              main:
                print foo
                print bar
              emit: 
                foo+bar
            }
            
            workflow delta() {
                take: foo
                take: bar
                main:
                println foo+bar
            }

            workflow empty { }
        '''

        when:
        def script = (TestScript)new GroovyShell(new ScriptBinding(), config).parse(SCRIPT).run()
        def meta = ScriptMeta.get(script)
        println meta.getWorkflow('bravo') .source.stripIndent()
        then:
        meta.definitions.size() == 4
        meta.getWorkflow('alpha') .declaredInputs == []
        meta.getWorkflow('alpha') .declaredVariables == []
        meta.getWorkflow('alpha') .source.stripIndent() == "print 'Hello world'\n"

        meta.getWorkflow('bravo') .declaredInputs == ['foo', 'bar']
        meta.getWorkflow('bravo') .declaredVariables == ['$out0']
        meta.getWorkflow('bravo') .source.stripIndent() == '''\
              take: foo
              take: bar
              main:
                print foo
                print bar
              emit: 
                foo+bar
              '''.stripIndent()

        meta.getWorkflow('delta') .declaredInputs == ['foo','bar']
        meta.getWorkflow('delta') .declaredVariables == [] 
        meta.getWorkflow('delta') .source.stripIndent() == '''\
                take: foo
                take: bar
                main:
                println foo+bar
                '''.stripIndent()

        meta.getWorkflow('empty') .source == ''
        meta.getWorkflow('empty') .declaredInputs == []
        meta.getWorkflow('empty') .declaredVariables == [] 
    }

    def 'should define anonymous workflow' () {

        def SCRIPT = '''
                    
            workflow {
              print 1
              print 2
            }
        '''

        when:
        def runner = new MockScriptRunner().setScript(SCRIPT).invoke()
        def meta = ScriptMeta.get(runner.getScript())
        then:
        meta.getWorkflow(null).getSource().stripIndent() == 'print 1\nprint 2\n'

    }

    def 'should run workflow block' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              take: foo
              emit: bar
              emit: baz  
              
              main: "$x world"
            }
       
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (TestScript)new GroovyShell(binding,config).parse(SCRIPT).run()
        def workflow = ScriptMeta.get(script).getWorkflow('alpha')
        then:
        workflow.declaredInputs == ['foo']
        workflow.declaredOutputs == ['bar', 'baz']

    }


    def 'should capture publish defs' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow {
              publish: 
                foo
                bar to: 'some/path'
                baz.out to: 'other/path'
              main:   
                x = 1 
            }
        '''

        when:
        def script = (TestScript)new GroovyShell(new ScriptBinding(), config).parse(SCRIPT).run()
        then:
        thrown(MultipleCompilationErrorsException)
    }

    def 'should not allow publish is sub-workflow' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              publish: foo
              main:   
                x = 1 
            }
        '''

        when:
        new GroovyShell(config).parse(SCRIPT)
        then:
        thrown(MultipleCompilationErrorsException)
    }

    def 'should report malformed workflow block' () {

        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              take: foo
              main: println foo
              take: bar
            }
       
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (TestScript)new GroovyShell(binding,config).parse(SCRIPT).run()
        def workflow = ScriptMeta.get(script).getWorkflow('alpha')
        then:
        def e = thrown(MultipleCompilationErrorsException)
        e.message.contains('Unexpected workflow `take` context here')

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
        def ch1 = new DataflowQueue(); ch1 << 'blah blah'
        def ch2 = new DataflowQueue(); ch2 << 'xxx'

        def binding = new WorkflowBinding(foo: 'Hello', bar: 'world', ch1: ch1, ch2: new ChannelOut([ch2]))
        def workflow = new WorkflowDef(binding: binding)

        when:
        def result = workflow.collectOutputs(['foo'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowVariable
        result.foo instanceof DataflowVariable
        result.foo.val == 'Hello'

        when:
        result = workflow.collectOutputs(['foo', 'bar'])
        then:
        result instanceof ChannelOut
        result.size()==2
        result[0] instanceof DataflowVariable
        result[1] instanceof DataflowVariable
        result.foo.val == 'Hello'
        result.bar.val == 'world'

        when:
        result = workflow.collectOutputs(['ch1'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowQueue
        result[0].val == 'blah blah'


        when:
        result = workflow.collectOutputs(['ch2'])
        then:
        result instanceof ChannelOut
        result.size()==1
        result[0] instanceof DataflowQueue
        result[0].val == 'xxx'

    }

    def 'should clone with a new name' () {
        given:
        def work = new WorkflowDef(name:'woo', body: new BodyDef({}, 'source'))

        when:
        def copy = work.cloneWithName('bar')
        then:
        copy.getName() == 'bar'
    }


    def 'should capture workflow code' () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
            workflow alpha {
              take:
                foo
              main:
                print x 
              emit: 
                foo  
            }
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (TestScript)new GroovyShell(binding,config).parse(SCRIPT).run()
        def workflow = ScriptMeta.get(script).getWorkflow('alpha')
        then:
        workflow.getSource().stripIndent() == '''\
                            take:
                              foo
                            main:
                              print x 
                            emit: 
                              foo  
                            '''.stripIndent()
    }

    def 'should capture empty workflow code'  () {
        given:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(TestScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
            workflow foo { } 
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (TestScript)new GroovyShell(binding,config).parse(SCRIPT).run()
        def workflow = ScriptMeta.get(script).getWorkflow('foo')
        then:
        workflow.getSource() == ''
    }

}
