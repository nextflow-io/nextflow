package nextflow.plugin

import com.sun.net.httpserver.HttpServer
import nextflow.Channel
import nextflow.extension.ChannelExtensionProvider
import test.Dsl2Spec
import test.MockScriptRunner

import java.nio.file.Files
import java.nio.file.Path

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginExtensionMethodsTest extends Dsl2Spec implements BuildPluginTrait{

    def folder

    def setup() {
        folder = buildPlugin( 'nf-plugin-template', '0.0.0', Path.of('../nf-plugin-template/build').toAbsolutePath())
    }

    def cleanup() {
        folder?.deleteDir()
        ChannelExtensionProvider.reset()
        Plugins.stop()
    }

    def 'should execute custom operator extension' () {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP

        where:
        SCRIPT_TEXT << ['''
            include { goodbye } from 'plugin/nf-plugin-template'

            channel
              .of('Bye bye folks')
              .goodbye()            
            ''', '''

            include { 
                reverse;
                goodbye; 
            } from 'plugin/nf-plugin-template'

            channel
              .of('Bye bye folks')
              .goodbye()             
            '''
        ]
    }

    def 'should execute custom factory extension' () {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

        where:
        SCRIPT_TEXT << ['''
            include { reverse } from 'plugin/nf-plugin-template'                

            channel.reverse('a string')            
            ''','''

            include { reverse;  } from 'plugin/nf-plugin-template'
            include { goodbye } from 'plugin/nf-plugin-template'
                
            channel.reverse('a string')            
            '''
        ]
    }

    def 'should execute custom operator as alias extension' () {
        given:
        def SCRIPT_TEXT = '''
            include { goodbye as myFunction } from 'plugin/nf-plugin-template'

            channel
              .of(100,200,300)
              .myFunction()            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP
    }

    def 'should execute custom factory as alias extension' () {
        given:
        def SCRIPT_TEXT = '''
            nextflow.enable.dsl=2
            include { reverse as myFunction } from 'plugin/nf-plugin-template'
         
            channel.myFunction('reverse this string')            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result
        result.val == 'reverse this string'.reverse()
        result.val == Channel.STOP

    }

    def 'should not include operators without the right signature' () {
        given:
        def SCRIPT_TEXT = '''
            nextflow.enable.dsl=2
            include { goodbyeWrongSignature } from 'plugin/nf-plugin-template'

            channel
              .of('Bye bye folks') | goodbyeWrongSignature                        
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(MissingMethodException)

    }

    def 'should not include factories without the right signature' () {
        given:
        def SCRIPT_TEXT = '''
            nextflow.enable.dsl=2
            include { reverseCantBeImportedBecauseWrongSignature } from 'plugin/nf-plugin-template'                

            channel.reverseCantBeImportedBecauseWrongSignature('a string')                        
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(IllegalStateException)

    }

}
