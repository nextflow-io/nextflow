package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Paths

import com.sun.net.httpserver.HttpServer
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginsFacadeTest extends Specification {

    def 'should setup plugins' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        def plugins = new PluginsFacade(folder)
        plugins.indexUrl = 'http://localhost:9900/plugins.json'

        when:
        plugins.setup([plugins: [ 'nf-console@1.0.0' ]])
        then:
        folder.resolve('nf-console-1.0.0').exists()

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
    }


    def 'should parse plugins config' () {
        given:
        def handler = new PluginsFacade()
        and:
        def cfg = [plugins: [ 'foo@1.2.3', 'bar@3.2.1' ]]

        when:
        final plugins = handler.parseConf(cfg)
        then:
        plugins.size() == 2
        and:
        plugins[0].id == 'foo'
        plugins[0].version == '1.2.3'
        and:
        plugins[1].id == 'bar'
        plugins[1].version == '3.2.1'
    }

    def 'should return plugin requirements' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-ignite': new PluginSpec('nf-ignite', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults, env: [:])

        when:
        def result = handler.pluginsRequirement([:])
        then:
        result == []

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT:'true'])
        result = handler.pluginsRequirement([:])
        then:
        result == [ new PluginSpec('nf-amazon', '0.1.0')]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT:'true'])
        result = handler.pluginsRequirement([tower:[enabled:true]])
        then:
        result == [
                new PluginSpec('nf-amazon', '0.1.0'),
                new PluginSpec('nf-tower', '0.1.0') ]
    }

    def 'should return default plugins given config' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-ignite': new PluginSpec('nf-ignite', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults)

        when:
        def plugins = handler.defaultPluginsConf([process:[executor: 'awsbatch']])
        then:
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-ignite' }
        
        when:
        plugins = handler.defaultPluginsConf([process:[executor: 'google-lifesciences']])
        then:
        plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-ignite' }

        when:
        plugins = handler.defaultPluginsConf([process:[executor: 'ignite']])
        then:
        plugins.find { it.id == 'nf-ignite' }
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }

        when:
        plugins = handler.defaultPluginsConf([:])
        then:
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-ignite' }
        !plugins.find { it.id == 'nf-google' }

    }

    def 'should get plugins list from env' () {

        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-ignite': new PluginSpec('nf-ignite', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT: 'nf-amazon,nf-tower@1.0.1,nf-foo@2.2.0,nf-bar'])

        when:
        def plugins = handler.defaultPluginsConf([:])
        then:
        plugins.size()==4
        plugins.find { it.id == 'nf-amazon' && it.version=='0.1.0' }    // <-- version from default
        plugins.find { it.id == 'nf-tower' && it.version=='1.0.1' }     // <-- version from the env var
        plugins.find { it.id == 'nf-foo' && it.version=='2.2.0' }       // <-- version from tne env var
        plugins.find { it.id == 'nf-bar' && it.version==null }          // <-- no version 
    }

    @Unroll
    def 'should validate plugins mode' () {
        given:
        def facade = new PluginsFacade(env: ENV)
        expect:
        facade.getPluginsMode() == EXPECTED
        where:
        ENV                             | EXPECTED
        [NXF_PLUGINS_MODE: 'foo']       | 'foo'
        [NXF_HOME: 'something']         | 'prod'
        [:]                             | 'dev'
    }

    @Unroll
    def 'should validate plugins default' () {
        given:
        def facade = new PluginsFacade(env: ENV)
        expect:
        facade.getPluginsDefault() == EXPECTED
        where:
        ENV                             | EXPECTED
        [NXF_PLUGINS_DEFAULT: 'true']   | true
        [NXF_PLUGINS_DEFAULT: 'nf-amzn']| true
        [NXF_PLUGINS_DEFAULT: 'false']  | false
        [NXF_HOME: 'something']         | true
        [:]                             | false
    }

    @Unroll
    def 'should validate plugins dir' () {
        given:
        def facade = new PluginsFacade(env: ENV)
        expect:
        facade.getPluginsDir() == EXPECTED
        where:
        ENV                                 | EXPECTED
        [NXF_PLUGINS_DIR: '/some/dir']      | Paths.get('/some/dir')
        [NXF_HOME: '/my/home']              | Paths.get('/my/home/plugins')
        [:]                                 | Paths.get('plugins')
    }
}
