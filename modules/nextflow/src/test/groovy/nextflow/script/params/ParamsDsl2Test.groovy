package nextflow.script.params

import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import spock.lang.Timeout
import test.Dsl2Spec
import test.MockScriptRunner
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ParamsDsl2Test extends Dsl2Spec {

    def 'should not allow unqualified input file' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: 
          tuple 'x'
          /touch x/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified input file declaration has been deprecated - replace `tuple 'x',..` with `tuple path('x'),..`"
    }

    def 'should not allow unqualified input val' () {
        given:
        def SCRIPT = '''
         
        process foo {
          input: 
          tuple X
          /echo $X/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified input value declaration has been deprecated - replace `tuple X,..` with `tuple val(X),..`"
    }


    def 'should not allow unqualified output file' () {
        given:
        def SCRIPT = '''
         
        process foo {
          output: 
          tuple 'x'
          /touch x/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified output path declaration has been deprecated - replace `tuple 'x',..` with `tuple path('x'),..`"
    }

    def 'should not allow unqualified output value' () {
        given:
        def SCRIPT = '''
         
        process foo {
          output: 
          tuple X
          /echo hello/
        }
       
        workflow {
            foo()
        }
        '''

        when:
        new MockScriptRunner() .setScript(SCRIPT).execute()
        then:
        def e = thrown(DeprecationException)
        e.message == "Unqualified output value declaration has been deprecated - replace `tuple X,..` with `tuple val(X),..`"
    }


    def 'should allow unqualified stdin and stdout' () {

        given:
        new Session()
        and:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
         process alpha {
              input:
              stdin
              output:
              stdout 

              /echo foo/
          }

         workflow { true } 
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (BaseScript)new GroovyShell(binding,config).parse(SCRIPT); script.run()
        and:
        def process = ScriptMeta.get(script).getProcess('alpha'); process.initialize()

        then:
        def inputs = process.processConfig.getInputs()
        def outputs = process.processConfig.getOutputs()
        and:
        inputs.size() == 1
        inputs[0] instanceof StdInParam
        and:
        outputs.size() == 1
        outputs[0] instanceof StdOutParam

    }

    def 'should allow unqualified tuple stdin and stdout' () {

        given:
        new Session()
        and:
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
                    
         process beta {
              input:
              tuple stdin, val(x)
              output:
              tuple stdout, path('z')

              /echo foo/
          }
       
         workflow { true } 
        '''

        when:
        def binding = new ScriptBinding().setSession(Mock(Session))
        def script = (BaseScript)new GroovyShell(binding,config).parse(SCRIPT); script.run()
        and:
        def process = ScriptMeta.get(script).getProcess('beta'); process.initialize()

        then:
        def inputs = process.processConfig.getInputs()
        def outputs = process.processConfig.getOutputs()
        and:
        inputs.size() == 1
        and:
        def in_tuple = (TupleInParam) inputs.get(0)
        and:
        in_tuple.inner.size() == 2
        in_tuple.inner[0] instanceof StdInParam
        in_tuple.inner[1] instanceof ValueInParam

        and:
        def out_tuple = (TupleOutParam) outputs.get(0)
        and:
        out_tuple.inner.size() == 2
        out_tuple.inner[0] instanceof StdOutParam
        out_tuple.inner[1] instanceof FileOutParam

    }
}
