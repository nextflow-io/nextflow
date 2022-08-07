package nextflow.plugin

import java.nio.file.Path

import nextflow.Channel
import nextflow.extension.ChannelExtensionProvider
import spock.lang.Shared
import test.Dsl2Spec
/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginExtensionMethodsTest extends Dsl2Spec {

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

    def 'should execute custom operator extension/1' () {
        given:
        def  SCRIPT_TEXT = '''
            include { goodbye } from 'plugin/nf-plugin-template'

            channel
              .of('Bye bye folks')
              .goodbye()            
            '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP
    }

    def 'should execute custom operator extension/2' () {
        given:
        def SCRIPT_TEXT = '''

            include { 
                reverse;
                goodbye; 
            } from 'plugin/nf-plugin-template'

            channel
              .of('Bye bye folks')
              .goodbye()             
            '''
        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP

    }

    def 'should execute custom factory extension/1' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverse } from 'plugin/nf-plugin-template'                

            channel.reverse('a string')            
            '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

    }

    def 'should execute custom factory extension/2' () {
        given:
        def SCRIPT_TEXT = '''

            include { reverse } from 'plugin/nf-plugin-template'
            include { goodbye } from 'plugin/nf-plugin-template'
                
            channel.reverse('a string')            
            '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

    }

    def 'should execute custom operator as alias extension' () {
        given:
        def SCRIPT_TEXT = '''
            include { goodbye as myFunction } from 'plugin/nf-plugin-template'

            channel
              .of(100,200,300)
              .myFunction()            
            '''

        when:
        def result = dsl_eval(SCRIPT_TEXT)

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

        when:
        def result = dsl_eval(SCRIPT_TEXT)

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

        when:
        dsl_eval(SCRIPT_TEXT)

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

        when:
        def result = dsl_eval(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)

    }

}
