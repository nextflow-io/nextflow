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
        plugins.env = [:]
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
        def defaults = new DefaultPlugins(plugins: [
                'delta': new PluginSpec('delta', '0.1.0'),
        ])

        def handler = new PluginsFacade(defaultPlugins: defaults)
        and:
        def cfg = [plugins: [ 'foo@1.2.3', 'bar@3.2.1', 'delta', 'omega' ]]

        when:
        final plugins = handler.parseConf(cfg)
        then:
        plugins.size() == 4
        and:
        plugins[0].id == 'foo'
        plugins[0].version == '1.2.3'
        and:
        plugins[1].id == 'bar'
        plugins[1].version == '3.2.1'
        and:
        plugins[2].id == 'delta'
        plugins[2].version == '0.1.0'
        and:
        plugins[3].id == 'omega'
        plugins[3].version == null
    }

    def 'should return plugin requirements' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
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
        result == []

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT:'true'])
        result = handler.pluginsRequirement([tower:[enabled:true]])
        then:
        result == [ new PluginSpec('nf-tower', '0.1.0') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [TOWER_ACCESS_TOKEN:'xyz'])
        result = handler.pluginsRequirement([:])
        then:
        result == [ new PluginSpec('nf-tower', '0.1.0') ]
    }

    def 'should return default plugins given config' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-azure': new PluginSpec('nf-azure', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults)

        when:
        def plugins = handler.defaultPluginsConf([process:[executor: 'awsbatch']])
        then:
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([process:[executor: 'google-lifesciences']])
        then:
        plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([process:[executor: 'azurebatch']])
        then:
        plugins.find { it.id == 'nf-azure' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([:])
        then:
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

    }

    def 'should return default plugins given workdir' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-azure': new PluginSpec('nf-azure', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults)

        when:
        def plugins = handler.defaultPluginsConf([workDir: 's3://foo'])
        then:
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([workDir: 'gs://foo'])
        then:
        plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([workDir: 'az://foo'])
        then:
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([:])
        then:
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

    }

    def 'should return default plugins given bucket dir' () {
        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-azure': new PluginSpec('nf-azure', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0')
        ])
        and:
        def handler = new PluginsFacade(defaultPlugins: defaults)

        when:
        def plugins = handler.defaultPluginsConf([bucketDir: 's3://foo'])
        then:
        plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([bucketDir: 'gs://foo'])
        then:
        plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([bucketDir: 'az://foo'])
        then:
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-amazon' }
        plugins.find { it.id == 'nf-azure' }

        when:
        plugins = handler.defaultPluginsConf([:])
        then:
        !plugins.find { it.id == 'nf-amazon' }
        !plugins.find { it.id == 'nf-google' }
        !plugins.find { it.id == 'nf-azure' }

    }

    def 'should get plugins list from env' () {

        given:
        def defaults = new DefaultPlugins(plugins: [
                'nf-amazon': new PluginSpec('nf-amazon', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
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


    // -- check Priority annotation

    @Priority(100)
    static class Foo { }

    def 'should get priority from object' () {
        given:
        def facade = new PluginsFacade()
        expect:
        facade.priority0('foo') == 0
        facade.priority0(new Foo()) == 100
    }

    // -- check priority ordering

    interface Bar { }

    @Priority(value = 10, group = 'alpha')
    static class XXX implements Bar {}

    @Priority(value = 20, group = 'alpha')
    static class YYY implements Bar {}

    @Priority(value = 30, group = 'beta')
    static class WWW implements Bar {}

    @Priority(40)
    static class ZZZ implements Bar {}

    def 'should get priority extensions' () {
        given:
        def xxx = new XXX()
        def yyy = new YYY()
        def zzz = new ZZZ()
        def www = new WWW()
        and:
        def THE_LIST = [yyy, xxx, zzz, www]; THE_LIST.shuffle()
        and:
        def facade = Spy( new PluginsFacade() ) {
            getExtensions(Foo) >> THE_LIST
        }

        when:
        def result = facade.getPriorityExtensions(Foo)
        then:
        // items are returned ordered by priority
        result.head() == xxx
        result.first() == xxx
        and:
        result == [xxx,yyy,www,zzz]

        when:
        result = facade.getPriorityExtensions(Foo, 'alpha')
        then:
        result == [xxx,yyy]

        when:
        result = facade.getPriorityExtensions(Foo, 'beta')
        then:
        result == [www]

    }

    // -- check scoped exceptions

    @Scoped(priority = 10, value = 'alpha')
    static class EXT1 implements Bar {}

    @Scoped(priority  = 20, value = 'alpha')
    static class EXT2 implements Bar {}

    @Scoped(priority  = 30, value = 'beta')
    static class EXT3 implements Bar {}

    @Scoped(priority  = -1, value = 'beta')
    static class EXT4 implements Bar {}

    @Scoped('')
    static class EXT5 implements Bar {}

    def 'should get scoped extensions' () {
        given:
        def ext1 = new EXT1()
        def ext2 = new EXT2()
        def ext3 = new EXT3()
        def ext4 = new EXT4()
        def ext5 = new EXT5()
        and:
        def THE_LIST = [ext1, ext2, ext3, ext4, ext5]; THE_LIST.shuffle()
        and:
        def facade = Spy( new PluginsFacade() ) {
            getExtensions(Foo) >> THE_LIST
        }

        when:
        def result = facade.getScopedExtensions(Foo)
        then:
        // items are returned ordered by priority
        result == [ ext1, ext4 ] as Set

        when:
        result = facade.getScopedExtensions(Foo, 'alpha')
        then:
        result == [ ext1 ] as Set

        when:
        result = facade.getScopedExtensions(Foo, 'beta')
        then:
        result == [ ext4 ] as Set

    }
}
