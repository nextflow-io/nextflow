package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Paths

import com.sun.net.httpserver.HttpServer
import nextflow.SysEnv
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
        plugins.init()
        plugins.load([plugins: [ 'nf-console@1.0.0' ]])
        then:
        folder.resolve('nf-console-1.0.0').exists()

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
    }

    def 'should create plugin manager' () {
        given:
        def folder = Files.createTempDirectory('test')
        def plugins = new PluginsFacade(folder,MODE)
        expect:
        plugins.createManager(folder,EMBEDDED).class == EXPECTED

        cleanup:
        folder?.deleteDir()

        where:
        MODE    | EMBEDDED      | EXPECTED
        'dev'   | false         | DevPluginManager
        'dev'   | true          | DevPluginManager
        and:
        'prod'  | false         | LocalPluginManager
        'prod'  | true          | EmbeddedPluginManager
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
                'nf-cloudcache': new PluginSpec('nf-cloudcache', '0.1.0'),
                'nf-google': new PluginSpec('nf-google', '0.1.0'),
                'nf-tower': new PluginSpec('nf-tower', '0.1.0'),
                'nf-wave': new PluginSpec('nf-wave', '0.1.0')
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
        result = handler.pluginsRequirement([tower:[enabled:false]])
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

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [:])
        result = handler.pluginsRequirement([wave:[enabled:true]])
        then:
        result == [ new PluginSpec('nf-wave', '0.1.0') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [:])
        result = handler.pluginsRequirement([plugins: [ 'foo@1.2.3']])
        then:
        result == [ new PluginSpec('foo', '1.2.3') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [:])
        result = handler.pluginsRequirement([plugins: [ 'nf-amazon@1.2.3']])
        then:
        result == [ new PluginSpec('nf-amazon', '1.2.3') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [:])
        result = handler.pluginsRequirement([plugins: [ 'nf-amazon']])
        then:
        result == [ new PluginSpec('nf-amazon', '0.1.0') ] // <-- config is taken from the default config

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT:'nf-google@2.0.0'])
        result = handler.pluginsRequirement([plugins: [ 'nf-amazon@1.2.3']])
        then:
        result == [ new PluginSpec('nf-amazon', '1.2.3'), new PluginSpec('nf-google','2.0.0') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [NXF_PLUGINS_DEFAULT:'nf-google@2.0.0'])
        result = handler.pluginsRequirement([:])
        then:
        result == [ new PluginSpec('nf-google','2.0.0') ]

        when:
        handler = new PluginsFacade(defaultPlugins: defaults, env: [:])
        result = handler.pluginsRequirement([cloudcache:[enabled:true]])
        then:
        result == [ new PluginSpec('nf-cloudcache', '0.1.0') ]

    }

    def 'should return default plugins given config' () {
        given:
        SysEnv.push([:])
        and:
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

        cleanup:
        SysEnv.pop()
    }

    def 'should return default plugins given workdir' () {
        given:
        SysEnv.push([:])
        and:
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

        cleanup:
        SysEnv.pop()
    }

    def 'should return default plugins given bucket dir' () {
        given:
        SysEnv.push([:])
        and:
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

        cleanup:
        SysEnv.pop()
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
        plugins.find { it.id == 'nf-foo' && it.version=='2.2.0' }       // <-- version from the env var
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

    @Unroll
    def 'should merge plugins' () {
        given:
        def facade = new PluginsFacade()
        def configPlugins = CONFIG.tokenize(',').collect { PluginSpec.parse(it) }
        def defaultPlugins = DEFAULT.tokenize(',').collect { PluginSpec.parse(it) }
        def expectedPlugins = EXPECTED.tokenize(',').collect { PluginSpec.parse(it) }

        expect:
        facade.mergePluginSpecs(configPlugins, defaultPlugins) == expectedPlugins

        where:
        CONFIG                  | DEFAULT               | EXPECTED
        ''                      | ''                    | ''
        'alpha,delta'           | ''                    | 'alpha,delta'
        ''                      | 'alpha,delta'         | 'alpha,delta'
        'alpha'                 | 'delta'               | 'alpha,delta'
        'delta'                 | 'alpha'               | 'delta,alpha'
        'delta'                 | 'delta'               | 'delta'
        'delta@1.0.0'           | 'delta'               | 'delta@1.0.0'
        'delta'                 | 'delta@2.0.0'         | 'delta@2.0.0'
        'delta@1.0.0'           | 'delta@2.0.0'         | 'delta@1.0.0'     // <-- config has priority
        'alpha,beta@1.0.0'      | 'delta@2.0.0'         | 'alpha,beta@1.0.0,delta@2.0.0'

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

}
