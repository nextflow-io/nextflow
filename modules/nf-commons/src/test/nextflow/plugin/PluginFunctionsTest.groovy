package nextflow.plugin

import com.sun.net.httpserver.HttpServer
import nextflow.Channel
import nextflow.NextflowMeta
import nextflow.exception.DuplicateModuleFunctionException
import nextflow.extension.ChannelExtensionProvider
import nextflow.extension.FunctionsExtensionProvider
import test.Dsl2Spec
import test.MockScriptRunner

import java.nio.file.Files

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginFunctionsTest extends Dsl2Spec {

    def 'should execute custom functions'() {
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
        FunctionsExtensionProvider.reset()
        Plugins.setup([plugins: ['nf-plugin-template@0.0.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == EXPECTED
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionProvider.reset()
        Plugins.stop()

        where:
        SCRIPT_TEXT                                                                           | EXPECTED
        "include { sayHello } from 'plugin/nf-plugin-template'; channel.of( sayHello() )"     | 'hi'
        "include { sayHello } from 'plugin/nf-plugin-template'; channel.of( sayHello('es') )" | 'hola'
        "include { sayHello as hi } from 'plugin/nf-plugin-template'; channel.of( hi() )"     | 'hi'

    }

    def 'should throw function not found'() {
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

        SCRIPT.text = '''
        include { sayHelloNotExist } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHelloNotExist() )
        '''

        when:
        FunctionsExtensionProvider.reset()
        Plugins.setup([plugins: ['nf-plugin-template@0.0.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(IllegalStateException)

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionProvider.reset()
        Plugins.stop()
    }

    def 'should not allow to include an existing function'() {
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

        SCRIPT.text = '''
        def sayHello(){ 'hi' }
        
        include { sayHello } from 'plugin/nf-plugin-template' 
        
        channel.of( sayHello() )
        '''

        when:
        FunctionsExtensionProvider.reset()
        Plugins.setup([plugins: ['nf-plugin-template@0.0.0']])

        NextflowMeta.instance.strictMode(true)
        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        thrown(DuplicateModuleFunctionException)

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionProvider.reset()
        Plugins.stop()
    }

    def 'should execute custom functions and channel extension at the same time'() {
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
        FunctionsExtensionProvider.reset()
        Plugins.setup([plugins: ['nf-plugin-template@0.0.0']])

        def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()

        then:
        result.val == EXPECTED
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionProvider.reset()
        Plugins.stop()

        where:
        SCRIPT_TEXT                                                                                           | EXPECTED
        "include { sayHello; goodbye } from 'plugin/nf-plugin-template'; channel.of( sayHello() ).goodbye() " | 'hi'

    }

}
