package nextflow.plugin

import nextflow.Channel
import nextflow.NextflowMeta
import nextflow.exception.DuplicateModuleFunctionException
import nextflow.extension.ChannelExtensionProvider

import test.Dsl2Spec
import test.MockScriptRunner

import java.nio.file.Files
import java.nio.file.Path


/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginFunctionsTest extends Dsl2Spec implements BuildPluginTrait{

    def folder

    def setup() {
        folder = buildPlugin( 'nf-plugin-template', '0.0.0', Path.of('../nf-plugin-template/build').toAbsolutePath())
    }

    def cleanup() {
        folder?.deleteDir()
        ChannelExtensionProvider.reset()
        Plugins.stop()
    }


    def 'should execute custom functions'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

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
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = '''
        include { sayHelloNotExist } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHelloNotExist() )
        '''

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(IllegalStateException)
    }

    def 'should not allow to include an existing function'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = '''
        def sayHello(){ 'hi' }
        
        include { sayHello } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHello() )
        '''

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(DuplicateModuleFunctionException)
    }

    def 'should allows to include an existing function but as alias'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = '''
        def sayHello(){ 'hi' }
        
        include { sayHello as anotherHello } from 'plugin/nf-plugin-template' 
        
        channel.of( anotherHello() )
        '''

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 'hi'
        result.val == Channel.STOP
    }

    def 'should not include a non annotated function'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = '''      
        include { aNonImportedFunction } from 'plugin/nf-plugin-template' 
        
        channel.of( aNonImportedFunction() )
        '''

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(IllegalStateException)
    }

    def 'should execute a function in a module'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
        nextflow.enable.dsl=2
                
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
        nextflow.enable.dsl=2
                
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
        nextflow.enable.dsl=2
                
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
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

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
        nextflow.enable.dsl=2
                
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
        nextflow.enable.dsl=2

        process sayHello {
            input:
              val lng
            output:
              stdout
              
            "$lng"
        }
        '''
        SCRIPT.text = """
        nextflow.enable.dsl=2
    
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
