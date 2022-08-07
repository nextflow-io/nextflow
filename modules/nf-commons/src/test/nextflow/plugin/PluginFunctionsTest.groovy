package nextflow.plugin


import java.nio.file.Path

import nextflow.Channel
import nextflow.exception.DuplicateModuleFunctionException
import nextflow.extension.ChannelExtensionProvider
import spock.lang.Shared
import spock.lang.TempDir
import test.Dsl2Spec
import test.MockScriptRunner
/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginFunctionsTest extends Dsl2Spec {

    @TempDir
    @Shared
    Path folder

    @Shared String pluginsMode

    def setup() {
        // this need to be set *before* the plugin manager class is created
        pluginsMode = System.getProperty('pf4j.mode')
        System.setProperty('pf4j.mode', 'dev')
        // the plugin root should
        def root = Path.of('.').toAbsolutePath().normalize()
        def manager = new TestPluginManager(root)
        Plugins.init(root, 'dev', manager)
    }

    def cleanup() {
        Plugins.stop()
        ChannelExtensionProvider.reset()
        pluginsMode ? System.setProperty('pf4j.mode',pluginsMode) : System.clearProperty('pf4j.mode')
    }

    def 'should execute custom functions'() {
        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result.val == EXPECTED
        result.val == Channel.STOP

        where:
        SCRIPT_TEXT                                                                           | EXPECTED
        "include { sayHello } from 'plugin/nf-plugin-template'; channel.of( sayHello() )"     | 'hi'
        "include { sayHello } from 'plugin/nf-plugin-template'; channel.of( sayHello('es') )" | 'hola'
        "include { sayHello as hi } from 'plugin/nf-plugin-template'; channel.of( hi() )"     | 'hi'

    }

    def 'should throw function not found'() {
        given:
        def SCRIPT_TEXT = '''
        include { sayHelloNotExist } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHelloNotExist() )
        '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)
    }

    def 'should not allow to include an existing function'() {
        given:
        def SCRIPT_TEXT = '''
        nextflow.enable.strict = true
        
        def sayHello(){ 'hi' }
        
        include { sayHello } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHello() )
        '''

        when:
        dsl_eval(SCRIPT_TEXT)

        then:
        thrown(DuplicateModuleFunctionException)
    }

    def 'should allows to include an existing function but as alias'() {
        given:
        def SCRIPT_TEXT= '''
        nextflow.enable.strict = true
        
        def sayHello(){ 'hi' }
        
        include { sayHello as anotherHello } from 'plugin/nf-plugin-template' 
        
        channel.of( anotherHello() )
        '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result.val == 'hi'
        result.val == Channel.STOP
    }

    def 'should not include a non annotated function'() {
        given:
        def SCRIPT_TEXT= '''      
        nextflow.enable.strict = true
        
        include { aNonImportedFunction } from 'plugin/nf-plugin-template' 
        
        channel.of( aNonImportedFunction() )
        '''

        when:
        dsl_eval(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)
    }

    def 'should execute a function in a module'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
        nextflow.enable.strict = true
                
        include { sayHello } from 'plugin/nf-plugin-template' 

        process foo {
            input:
              val lng
            output:
              stdout
              
            "${sayHello(lng)}"
        }
        '''

        SCRIPT.text = '''
        include { foo } from './module.nf'        
        workflow{
            main:
                foo( 'en' )
            emit:
                foo.out
        }
        '''

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 'hi'

    }

    def 'should execute a function in two modules'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE1 = folder.resolve('module1.nf')
        def MODULE2 = folder.resolve('module2.nf')

        MODULE1.text = '''
        nextflow.enable.strict=true
                
        include { sayHello } from 'plugin/nf-plugin-template' 

        process foo {
            input:
              val lng
            output:
              stdout
              
            "${sayHello(lng)}"
        }
        '''

        MODULE2.text = '''
        nextflow.enable.strict=true
                
        include { sayHello } from 'plugin/nf-plugin-template' 

        process bar {
            input:
              val lng
            output:
              stdout
              
            "${sayHello('es')}"
        }
        '''

        SCRIPT.text = '''
        include { foo } from './module1.nf'        
        include { bar } from './module2.nf'
        
        workflow{
            main:
                foo( 'en' ) | bar
            emit:
                bar.out
        }
        '''

        when:
        def result = dsl_eval(SCRIPT)

        then:
        result.val == 'hola'
    }

    def 'should execute custom functions and channel extension at the same time'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == EXPECTED
        result.val == Channel.STOP

        where:
        SCRIPT_TEXT                                                                                           | EXPECTED
        "include { sayHello; goodbye } from 'plugin/nf-plugin-template'; channel.of( sayHello() ).goodbye() " | 'hi'

    }

    def 'should not allow a function with the same name as a process'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = """
        nextflow.enable.strict=true
                
        include { sayHello } from 'plugin/nf-plugin-template' 

        process sayHello {
            input:
              val lng
            output:
              stdout
              
            "Hi"
        }
        workflow{
            main:
                Channel.from('hi') | sayHello
            emit:
                sayHello.out
        }
        """.stripIndent()

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(DuplicateModuleFunctionException)
    }

    def 'should not allow a function and a process with the same from other module'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE1 = folder.resolve('module1.nf')

        MODULE1.text = '''
        nextflow.enable.strict=true

        process sayHello {
            input:
              val lng
            output:
              stdout
              
            "$lng"
        }
        '''
        SCRIPT.text = """
        nextflow.enable.strict=true
    
        include { sayHello } from 'plugin/nf-plugin-template'                 
        include { sayHello } from './module1.nf' 

        process foo {
            input:
              val lng
            output:
              stdout
              
            "Hi"
        }
        workflow{
            main:
                Channel.from('hi') | foo
            emit:
                foo.out
        }
        """.stripIndent()

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(DuplicateModuleFunctionException)

    }

}
