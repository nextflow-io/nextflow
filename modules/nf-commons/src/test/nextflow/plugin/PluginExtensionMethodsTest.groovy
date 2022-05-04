package nextflow.plugin

import com.sun.net.httpserver.HttpServer
import nextflow.Channel
import nextflow.extension.ChannelExtensionProvider
import test.Dsl2Spec
import test.MockScriptRunner

import java.nio.file.Files

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginExtensionMethodsTest extends Dsl2Spec {

    def 'should execute custom operator extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()

        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-hello@0.2.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionProvider.reset()
        Plugins.stop()

        where:
        SCRIPT_TEXT << ['''
            include { goodbye } from 'plugin/nf-hello'

            channel
              .of('Bye bye folks')
              .goodbye()            
            ''', '''

            include { 
                reverse;
                goodbye; 
            } from 'plugin/nf-hello'

            channel
              .of('Bye bye folks')
              .goodbye()             
            '''
        ]
    }

    def 'should execute custom factory extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-hello@0.2.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionProvider.reset()

        where:
        SCRIPT_TEXT << ['''
            include { reverse } from 'plugin/nf-hello'                

            channel.reverse('a string')            
            ''','''

            include { reverse;  } from 'plugin/nf-hello'
            include { goodbye } from 'plugin/nf-hello'
                
            channel.reverse('a string')            
            '''
        ]
    }

    def 'should execute custom operator as alias extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def SCRIPT_TEXT = '''
            include { goodbye as myFunction } from 'plugin/nf-hello'

            channel
              .of(100,200,300)
              .myFunction()            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-hello@0.2.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionProvider.reset()
    }

    def 'should execute custom factory as alias extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def SCRIPT_TEXT = '''
            nextflow.enable.dsl=2
            include { reverse as myFunction } from 'plugin/nf-hello'
         
            channel.myFunction('reverse this string')            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-hello@0.2.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result
        result.val == 'reverse this string'.reverse()
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionProvider.reset()
    }

}
