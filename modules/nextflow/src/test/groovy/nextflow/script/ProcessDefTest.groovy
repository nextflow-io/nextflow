package nextflow.script

import spock.lang.Specification

import nextflow.Session
import nextflow.ast.NextflowDSL
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import test.MockSession
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDefTest extends Specification {

    def 'should define processes' () {

        given:
        def session = new MockSession()
        def binding = new ScriptBinding(session)
        def config = new CompilerConfiguration()
        config.setScriptBaseClass(BaseScript.class.name)
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        def SCRIPT = '''
            
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
            
        '''

        when:
        def script = (BaseScript)new GroovyShell(binding,config).parse(SCRIPT)

        then:
        true

    }


    def 'should clone a process with a new name '() {

        given:
        def OWNER = Mock(BaseScript)
        def BODY = { -> null }
        def proc = new ProcessDef(OWNER, BODY, 'foo')

        when:
        def copy = proc.cloneWithName('foo_alias')
        then:
        copy.getName() == 'foo_alias'
        copy.getSimpleName() == 'foo_alias'
        copy.getBaseName() == 'foo'
        copy.getOwner() == OWNER
        copy.rawBody.class == BODY.class
        !copy.rawBody.is(BODY)

        when:
        copy = proc.cloneWithName('flow1:flow2:foo')
        then:
        copy.getName() == 'flow1:flow2:foo'
        copy.getSimpleName() == 'foo'
        copy.getBaseName() == 'foo'
        copy.getOwner() == OWNER
        copy.rawBody.class == BODY.class
        !copy.rawBody.is(BODY)
    }

    def 'should apply process config' () {
        given:
        def OWNER = Mock(BaseScript)
        def CONFIG = [
                process:[
                        cpus:2, memory: '1GB',
                        'withName:foo': [memory: '3GB'],
                        'withName:bar': [cpus:4, memory: '4GB'],
                        'withName:flow1:flow2:flow3:bar': [memory: '8GB']
                ]
        ]
        def BODY = {->
            return new BodyDef({->}, 'echo hello')
        }
        def proc = new ProcessDef(OWNER, BODY, 'foo')
        proc.session = Mock(Session) { getConfig() >> CONFIG }

        when:
        def copy = proc.clone()
        copy.initialize()
        then:
        def cfg1 = copy.processConfig.createTaskConfig()
        cfg1.getCpus()==2           // taken from the generic config
        cfg1.getMemory().giga == 3  // taken from the `foo` config

        when:
        copy = proc.cloneWithName('flow1:bar')
        copy.initialize()
        then:
        def cfg2 = copy.processConfig.createTaskConfig()
        cfg2.getCpus()==4           // taken from the `bar` config
        cfg2.getMemory().giga == 4  // taken from the `bar` config


        when:
        copy = proc.cloneWithName('flow1:flow2:flow3:bar')
        copy.initialize()
        then:
        def cfg3 = copy.processConfig.createTaskConfig()
        cfg3.getCpus()==4           // <-- taken from `withName: foo`
        cfg3.getMemory().giga == 8  // <-- taken from `withName: 'flow1:flow2:flow3:bar'`
    }
}
