/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.plugin

import static com.github.tomakehurst.wiremock.client.WireMock.*

import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

import ch.qos.logback.classic.Logger
import ch.qos.logback.classic.spi.ILoggingEvent
import ch.qos.logback.core.read.ListAppender
import com.github.tomakehurst.wiremock.junit.WireMockRule
import nextflow.exception.AbortOperationException
import org.junit.Rule
import org.pf4j.ExtensionPoint
import org.pf4j.PluginRuntimeException
import org.pf4j.PluginState
import org.slf4j.LoggerFactory
import spock.lang.IgnoreIf
import spock.lang.Specification

/**
 * A moved plugin keeps `Plugin-Id: nf-foo` and declares `Plugin-Namespace: myorg`
 */
class PluginNamespaceTest extends Specification {

    @Rule
    WireMockRule wiremock = new WireMockRule(0)

    Path folder
    Path repo
    Path store
    int count

    static class FooExtension implements ExtensionPoint { }

    def setup() {
        // other tests leave pf4j in dev mode, where core plugins are not downloaded
        System.setProperty('pf4j.mode', 'deployment')
        folder = Files.createTempDirectory('test')
        repo = Files.createDirectory(folder.resolve('repo'))
        store = Files.createDirectory(folder.resolve('plugins'))
    }

    def cleanup() {
        store.toFile().setWritable(true)
        folder?.deleteDir()
    }

    // as the Gradle plugin builds it: `Plugin-Id` is the plain name
    private Path zip(Map opts=[:], String name, String ver) {
        final dir = repo.resolve("z${count++}").resolve("$name-$ver")
        final meta = Files.createDirectories(dir.resolve('classes/META-INF'))
        meta.resolve('MANIFEST.MF').text = """\
            Manifest-Version: 1.0
            Plugin-Class: ${PluginUpdaterTest.FooPlugin.name}
            Plugin-Id: ${opts.id ?: name}
            Plugin-Version: $ver
            """.stripIndent() + (opts.ns ? "Plugin-Namespace: $opts.ns\n" : '') + (opts.deps ? "Plugin-Dependencies: $opts.deps\n" : '')
        if( opts.ext )
            meta.resolve('extensions.idx').text = FooExtension.name + '\n'
        final result = dir.resolveSibling("${dir.fileName}.zip")
        new ZipOutputStream(Files.newOutputStream(result)).withCloseable { zos ->
            for( Path f : Files.list(meta).toList() ) {
                zos.putNextEntry(new ZipEntry(dir.relativize(f).toString()))
                zos.write(f.bytes)
                zos.closeEntry()
            }
        }
        return result
    }

    // an unzipped plugin, as stored by an older Nextflow
    private Path stored(Map opts=[:], String folderName, String name, String ver) {
        final target = store.resolve(folderName)
        final dir = zip(opts, name, ver).resolveSibling("$name-$ver")
        Files.move(dir, target)
        return target
    }

    private static String entry(String id, String canonicalId, Path... zips) {
        final releases = zips.collect { z -> """{"version":"${z.fileName.toString().replaceAll(/.*-(\d+\.\d+\.\d+)\.zip/, '$1')}","url":"file:$z","date":"2025-03-21T10:40:36Z"}""" }
        return """{"id":"$id",${canonicalId ? "\"canonicalId\":\"$canonicalId\"," : ''}"releases":[${releases.join(',')}]}"""
    }

    private void registry(String... entries) {
        wiremock.stubFor(get(urlPathEqualTo('/v1/plugins/dependencies'))
            .willReturn(aResponse().withStatus(200).withBody("""{"plugins":[${entries.join(',')}]}""")))
    }

    private int registryCalls() {
        wiremock.findAll(getRequestedFor(urlPathEqualTo('/v1/plugins/dependencies'))).size()
    }

    private URL getRegistryUrl() { new URL("http://localhost:${wiremock.port()}/") }

    private PluginUpdater updater(boolean offline=false) {
        new PluginUpdater(new LocalPluginManager(store), store, offline ? null : registryUrl, offline)
    }

    private PluginsFacade facade(Map env=[:]) {
        final result = new PluginsFacade(store, 'prod', false, registryUrl.toString())
        result.@env = env
        result.init()
        return result
    }

    private static List<String> loaded(PluginUpdater updater) {
        updater.@pluginManager.plugins*.pluginId.sort()
    }

    private boolean isLink(String name, String target) {
        final path = store.resolve(name)
        Files.isSymbolicLink(path) && Files.readSymbolicLink(path).toString() == target
    }

    private boolean isFolder(String name) {
        Files.isDirectory(store.resolve(name), LinkOption.NOFOLLOW_LINKS)
    }

    // -- loading

    def 'should load a moved plugin once, by both IDs' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: NS)
        registry(entry('nf-foo', 'myorg/nf-foo', z), entry('myorg/nf-foo', 'myorg/nf-foo', z))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('nf-foo', '1.0.0')])
        updater.prepareAndStart('nf-foo', '1.0.0')
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
        updater.@pluginManager.getPlugin('myorg/nf-foo').pluginState == PluginState.STARTED
        and: 'stored once, linked under the plain name'
        isFolder('myorg~nf-foo-1.0.0')
        isLink('nf-foo-1.0.0', 'myorg~nf-foo-1.0.0')

        where:
        NS << ['myorg', null]
    }

    def 'should learn the namespace from the zip of a plain request' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', null, z))
        def updater = updater()

        when:
        updater.prepareAndStart('nf-foo', '1.0.0')
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
        updater.aliases == ['nf-foo': 'myorg/nf-foo']
        isFolder('nf-foo-1.0.0')
        isLink('myorg~nf-foo-1.0.0', 'nf-foo-1.0.0')
    }

    def 'should reuse the loaded plugin when the plain name comes after the namespaced ID' () {
        given: 'the registry gives no canonicalId'
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', null, z), entry('myorg/nf-foo', null, z))
        def updater = updater()

        when:
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')
        updater.prepareAndStart('nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
        updater.aliases == ['nf-foo': 'myorg/nf-foo']
    }

    def 'should load a zip without Plugin-Namespace by the namespaced ID' () {
        given:
        def z = zip('nf-foo', '1.0.0')
        registry(entry('myorg/nf-foo', 'myorg/nf-foo', z))
        def updater = updater()

        when:
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then: 'the namespace comes from the store folder name'
        loaded(updater) == ['myorg/nf-foo']
        isFolder('myorg~nf-foo-1.0.0')
    }

    def 'should load a zip with a namespaced Plugin-Id' () {
        given:
        def z = zip('nf-foo', '1.0.0', id: 'myorg/nf-foo')
        registry(entry('myorg/nf-foo', null, z))
        def updater = updater()

        when:
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
    }

    def 'should refuse a zip that declares another namespace' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'other')
        registry(entry('myorg/nf-foo', null, z))
        def updater = updater()

        when:
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        def e = thrown(PluginRuntimeException)
        e.message == 'Plugin myorg/nf-foo cannot be loaded -- its zip declares the ID other/nf-foo'
        loaded(updater) == []
    }

    def 'should load the same name from two namespaces' () {
        given:
        def z1 = zip('nf-foo', '1.0.0', ns: 'org1')
        def z2 = zip('nf-foo', '1.0.0', ns: 'org2')
        registry(entry('org1/nf-foo', null, z1), entry('org2/nf-foo', null, z2))
        def updater = updater()

        when:
        updater.prepareAndStart('org1/nf-foo', '1.0.0')
        updater.prepareAndStart('org2/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['org1/nf-foo', 'org2/nf-foo']
    }

    def 'should not change an unscoped plugin that did not move' () {
        given:
        def z = zip('nf-bar', '1.0.0')
        registry(entry('nf-bar', null, z))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('nf-bar', '1.0.0')])
        updater.prepareAndStart('nf-bar', '1.0.0')

        then:
        loaded(updater) == ['nf-bar']
        updater.aliases.isEmpty()
        store.toFile().list() as List == ['nf-bar-1.0.0']
        registryCalls() == 1

        when: 'installed and pinned'
        this.updater().prefetchMetadata([new PluginRef('nf-bar', '1.0.0')])
        then:
        registryCalls() == 1
    }

    // -- canonicalId from the registry

    def 'should map the plugins named in the registry response' () {
        given:
        def bar = zip('nf-bar', '1.0.0', deps: 'nf-foo@>=1.0.0')
        def foo = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-bar', null, bar), entry('nf-foo', CANONICAL, foo))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('nf-bar', '1.0.0')])

        then:
        updater.aliases == EXPECTED
        updater.getPluginsMap().keySet() == KEYS as Set

        where:
        CANONICAL       | EXPECTED                      | KEYS
        'myorg/nf-foo'  | ['nf-foo': 'myorg/nf-foo']    | ['nf-bar', 'nf-foo', 'myorg/nf-foo']
        null            | [:]                           | ['nf-bar', 'nf-foo']
    }

    def 'should find a moved plugin by namespaced ID in a legacy plugins.json' () {
        given: 'a plugins.json index and a store link'
        def z = zip('nf-foo', '1.0.0')
        def index = repo.resolve('plugins.json')
        index.text = """[{"id":"nf-foo","releases":[{"version":"1.0.0","date":"Jun 25, 2020 9:58:35 PM","url":"file:$z"}]}]"""
        stored('nf-foo-0.9.0', 'nf-foo', '0.9.0')
        Files.createSymbolicLink(store.resolve('myorg~nf-foo-0.9.0'), Path.of('nf-foo-0.9.0'))
        def updater = new PluginUpdater(new LocalPluginManager(store), store, index.toUri().toURL(), false)

        when:
        updater.prepareAndStart('nf-foo', '1.0.0')

        then:
        updater.getPluginsMap().get('myorg/nf-foo').id == 'nf-foo'
        loaded(updater) == ['myorg/nf-foo']
        isLink('nf-foo-1.0.0', 'myorg~nf-foo-1.0.0')
    }

    def 'should map a plain dependency on a moved plugin' () {
        given:
        def bar = zip('nf-bar', '1.0.0', deps: 'nf-foo?@>=1.0.0')
        def foo = zip('nf-foo', '1.0.0')
        and: 'LATE: the canonicalId comes only with the dependency metadata, after nf-bar is loaded'
        wiremock.stubFor(get(urlPathEqualTo('/v1/plugins/dependencies')).withQueryParam('plugins', matching('nf-bar.*'))
            .willReturn(aResponse().withStatus(200).withBody("""{"plugins":[${entry('nf-bar', null, bar)}${LATE ? '' : ',' + entry('nf-foo', 'myorg/nf-foo', foo)}]}""")))
        wiremock.stubFor(get(urlPathEqualTo('/v1/plugins/dependencies')).withQueryParam('plugins', matching('nf-foo.*'))
            .willReturn(aResponse().withStatus(200).withBody("""{"plugins":[${entry('nf-foo', 'myorg/nf-foo', foo)}]}""")))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('nf-bar', '1.0.0')])
        updater.prepareAndStart('nf-bar', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo', 'nf-bar']
        updater.@pluginManager.getPlugin('nf-bar').pluginState == PluginState.STARTED
        with(updater.@pluginManager.getPlugin('nf-bar').descriptor.dependencies[0]) {
            pluginId == 'myorg/nf-foo'
            optional
            pluginVersionSupport == '>=1.0.0'
        }

        where:
        LATE << [false, true]
    }

    def 'should keep the pinned version of a moved core plugin' () {
        given:
        def bar = zip('nf-bar', '1.0.0', deps: 'nf-core@>=1.0.0')
        def core1 = zip('nf-core', '1.0.0', ns: 'nextflow-io')
        def core2 = zip('nf-core', '2.0.0', ns: 'nextflow-io')
        registry(entry('nf-bar', null, bar), entry('nf-core', 'nextflow-io/nf-core', core1, core2))
        def updater = updater()
        updater.@defaultPlugins = new DefaultPlugins(plugins: ['nextflow-io/nf-core': new PluginRef('nextflow-io/nf-core', '1.0.0')])

        when:
        updater.prefetchMetadata([new PluginRef('nf-bar', '1.0.0')])
        updater.prepareAndStart('nf-bar', '1.0.0')

        then:
        loaded(updater) == ['nextflow-io/nf-core', 'nf-bar']
        updater.@pluginManager.getPlugin('nextflow-io/nf-core').descriptor.version == '1.0.0'
    }

    // -- call sites

    def 'should start a moved plugin once when a run names both IDs' () {
        given:
        def z1 = zip('nf-foo', '1.0.0', ns: 'myorg')
        def z2 = zip('nf-foo', '2.0.0', ns: 'myorg')
        registry(entry('nf-foo', 'myorg/nf-foo', z1, z2), entry('myorg/nf-foo', 'myorg/nf-foo', z1, z2))
        def plugins = facade()

        when: 'the config ID comes first, as pluginsRequirement lists it'
        plugins.start([new PluginRef('myorg/nf-foo', '1.0.0'), new PluginRef('nf-foo', '2.0.0')])

        then:
        plugins.manager.plugins*.pluginId == ['myorg/nf-foo']
        plugins.getPlugin('nf-foo').descriptor.version == '1.0.0'
        !store.resolve('myorg~nf-foo-2.0.0').exists()
    }

    def 'should let the config version win over the default plugin version' () {
        given:
        def z1 = zip('nf-foo', '1.0.0', ns: 'myorg')
        def z2 = zip('nf-foo', '2.0.0', ns: 'myorg')
        registry(entry('nf-foo', 'myorg/nf-foo', z1, z2), entry('myorg/nf-foo', 'myorg/nf-foo', z1, z2))
        def plugins = facade(NXF_PLUGINS_DEFAULT: 'nf-foo@2.0.0')

        when:
        plugins.load([plugins: ['myorg/nf-foo@1.0.0']])

        then:
        plugins.manager.plugins*.pluginId == ['myorg/nf-foo']
        plugins.getPlugin('nf-foo').descriptor.version == '1.0.0'
    }

    def 'should resolve the plain ID of an include' () {
        given:
        def z = zip('nf-foo', '1.0.0')
        registry(entry('nf-foo', 'myorg/nf-foo', z), entry('myorg/nf-foo', 'myorg/nf-foo', z))

        when: 'the namespaced plugin is started'
        def plugins = facade()
        plugins.start('myorg/nf-foo@1.0.0')
        then:
        !plugins.startIfMissing('nf-foo')
        plugins.isStarted('nf-foo')
        plugins.manager.plugins*.pluginId == ['myorg/nf-foo']

        when: 'nothing is started'
        plugins = facade()
        then:
        plugins.startIfMissing('nf-foo')
        plugins.isStarted('nf-foo')
        plugins.isStarted('myorg/nf-foo')
        plugins.manager.plugins*.pluginId == ['myorg/nf-foo']
    }

    def 'should find a moved plugin and its extensions by the plain ID' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg', ext: true)
        registry(entry('myorg/nf-foo', 'myorg/nf-foo', z))
        def plugins = facade()

        when:
        plugins.start('myorg/nf-foo@1.0.0')

        then:
        plugins.getPlugin('nf-foo').pluginId == 'myorg/nf-foo'
        plugins.getExtensions(FooExtension, 'nf-foo').size() == 1
        and: 'an unknown plain name is not mapped'
        plugins.getPlugin('nf-bar') == null
    }

    def 'should check the allowed plugins by both IDs' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', 'myorg/nf-foo', z), entry('myorg/nf-foo', 'myorg/nf-foo', z))
        def plugins = facade(NXF_PLUGINS_ALLOWED: ALLOWED)

        when:
        plugins.start(REQUEST)

        then:
        plugins.getPlugin(REQUEST.tokenize('@')[0]).pluginId == 'myorg/nf-foo'

        where:
        ALLOWED         | REQUEST
        'myorg/nf-foo'  | 'nf-foo@1.0.0'
        'nf-foo'        | 'myorg/nf-foo@1.0.0'
    }

    def 'should refuse a plugin that is not allowed' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', null, z))
        def plugins = facade(NXF_PLUGINS_ALLOWED: 'myorg/nf-foo')

        when: 'no canonicalId says nf-foo moved'
        plugins.start('nf-foo@1.0.0')

        then:
        thrown(AbortOperationException)
        plugins.manager.plugins.isEmpty()
    }

    def 'should warn once when the config names a moved plugin' () {
        given:
        def logger = (Logger) LoggerFactory.getLogger(PluginsFacade)
        def appender = new ListAppender<ILoggingEvent>(); appender.start()
        logger.addAppender(appender)
        and:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', 'myorg/nf-foo', z), entry('myorg/nf-foo', 'myorg/nf-foo', z))
        def plugins = facade()

        when:
        plugins.load([plugins: ['nf-foo@1.0.0']])
        plugins.load([plugins: ['nf-foo@1.0.0']])

        then:
        appender.list*.formattedMessage.findAll { it.contains('has moved') } == ['Plugin `nf-foo` has moved to `myorg/nf-foo` -- update your config']
        plugins.manager.plugins*.pluginId == ['myorg/nf-foo']

        cleanup:
        logger.detachAppender(appender)
    }

    def 'should find an embedded plugin by its plain ID' () {
        given:
        stored('nf-foo-1.0.0', 'nf-foo', '1.0.0', ns: 'myorg')
        def plugins = new PluginsFacade(store, 'prod', false, registryUrl.toString())

        when:
        plugins.init(true)

        then:
        plugins.getPlugin('nf-foo').pluginId == 'myorg/nf-foo'
        plugins.isStarted('nf-foo')
    }

    // -- plugin store

    def 'should reuse the plain folder of an older Nextflow' () {
        given: 'the registry release cannot be downloaded'
        stored('nf-foo-1.0.0', 'nf-foo', '1.0.0')
        registry(entry('myorg/nf-foo', 'myorg/nf-foo', folder.resolve('missing/nf-foo-1.0.0.zip')))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('myorg/nf-foo', '1.0.0')])
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
        isLink('myorg~nf-foo-1.0.0', 'nf-foo-1.0.0')
    }

    def 'should never link to a link or to an incomplete folder' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'org1')
        registry(entry('org1/nf-foo', 'org1/nf-foo', z), entry('nf-foo', 'org1/nf-foo', z))
        and:
        if( LINK ) {
            stored('org2~nf-foo-1.0.0', 'nf-foo', '1.0.0', ns: 'org2')
            Files.createSymbolicLink(store.resolve('nf-foo-1.0.0'), Path.of('org2~nf-foo-1.0.0'))
        }
        else
            Files.createDirectory(store.resolve('nf-foo-1.0.0'))
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('org1/nf-foo', '1.0.0')])
        updater.prepareAndStart('org1/nf-foo', '1.0.0')

        then: 'downloaded, not linked'
        loaded(updater) == ['org1/nf-foo']
        isFolder('org1~nf-foo-1.0.0')

        where:
        LINK << [true, false]
    }

    @IgnoreIf({ System.getProperty('user.name') == 'root' })
    def 'should load the plain folder of an older Nextflow from a read-only store' () {
        given: 'the registry release cannot be downloaded'
        stored('nf-foo-1.0.0', 'nf-foo', '1.0.0')
        registry(entry('myorg/nf-foo', 'myorg/nf-foo', folder.resolve('missing/nf-foo-1.0.0.zip')))
        store.toFile().setWritable(false)
        def updater = updater()

        when:
        updater.prefetchMetadata([new PluginRef('myorg/nf-foo', '1.0.0')])
        updater.prepareAndStart('myorg/nf-foo', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo']
        !Files.exists(store.resolve('myorg~nf-foo-1.0.0'), LinkOption.NOFOLLOW_LINKS)
    }

    def 'should refuse an invalid Plugin-Namespace' () {
        given:
        def dir = stored('nf-foo-1.0.0', 'nf-foo', '1.0.0', ns: NS)

        when:
        new NamespaceDescriptorFinder().find(dir)

        then:
        def e = thrown(PluginRuntimeException)
        e.message == "Invalid Plugin-Namespace '$NS' in plugin $dir"

        where:
        NS << ['a~b', 'a/b', '..', '.']
    }

    def 'should map only from store links' () {
        given:
        ['org1~nf-foo-1.0.0', 'nf-bar-1.2.0', 'org1~nf-qux-1.0.0', 'org2~nf-qux-1.0.0', 'org1~nf-baz-1.0.0', 'nf-baz-1.0.0'].each { Files.createDirectory(store.resolve(it)) }
        Files.createSymbolicLink(store.resolve('nf-foo-1.0.0'), Path.of('org1~nf-foo-1.0.0'))
        Files.createSymbolicLink(store.resolve('org2~nf-bar-1.2.0'), Path.of('nf-bar-1.2.0'))
        and: 'a name in two namespaces, and a broken link'
        Files.createSymbolicLink(store.resolve('nf-qux-1.0.0'), Path.of('org1~nf-qux-1.0.0'))
        Files.createSymbolicLink(store.resolve('nf-qux-2.0.0'), Path.of('org2~nf-qux-1.0.0'))
        Files.createSymbolicLink(store.resolve('nf-gone-1.0.0'), Path.of('org1~nf-gone-1.0.0'))

        expect: 'nf-baz has folders with both names, but no link'
        PluginUpdater.storeAliases(store) == ['nf-foo': 'org1/nf-foo', 'nf-bar': 'org2/nf-bar']
    }

    def 'should ask the registry when a run names both IDs and the store has no link' () {
        given:
        def z = zip('nf-foo', '1.0.0', ns: 'myorg')
        registry(entry('nf-foo', 'myorg/nf-foo', z), entry('myorg/nf-foo', 'myorg/nf-foo', z))
        stored('nf-foo-1.0.0', 'nf-foo', '1.0.0', ns: 'myorg')
        stored('myorg~nf-foo-1.0.0', 'nf-foo', '1.0.0', ns: 'myorg')
        def both = [new PluginRef('nf-foo', '1.0.0'), new PluginRef('myorg/nf-foo', '1.0.0')]

        when:
        updater().prefetchMetadata(both)
        then:
        registryCalls() == 1

        when: 'one name alone'
        updater().prefetchMetadata(both.take(1))
        then:
        registryCalls() == 1

        when: 'a link says they are one plugin'
        store.resolve('nf-foo-1.0.0').deleteDir()
        Files.createSymbolicLink(store.resolve('nf-foo-1.0.0'), Path.of('myorg~nf-foo-1.0.0'))
        updater().prefetchMetadata(both)
        then:
        registryCalls() == 1
    }

    def 'should resolve offline a plain dependency on a moved plugin' () {
        given: 'a store filled online, with the canonicalId'
        def foo = zip('nf-foo', '1.0.0')
        def bar = zip('nf-bar', '1.0.0', deps: 'nf-foo@1.0.0')
        registry(entry('nf-bar', null, bar), entry('nf-foo', 'myorg/nf-foo', foo), entry('myorg/nf-foo', 'myorg/nf-foo', foo))
        updater().pullPlugins(['nf-bar@1.0.0'])
        assert isFolder('myorg~nf-foo-1.0.0') && isLink('nf-foo-1.0.0', 'myorg~nf-foo-1.0.0')
        and: 'a broken link'
        Files.createSymbolicLink(store.resolve('nf-gone-1.0.0'), Path.of('myorg~nf-gone-1.0.0'))
        def updater = updater(true)

        when:
        updater.prepareAndStart('nf-bar', '1.0.0')

        then:
        loaded(updater) == ['myorg/nf-foo', 'nf-bar']
        updater.@pluginManager.getPlugin('nf-bar').pluginState == PluginState.STARTED
    }

    def 'should add the missing links on plugin install' () {
        given: 'a store filled before the move'
        def z = zip('nf-foo', '1.0.0')
        registry(entry('myorg/nf-foo', null, z))
        updater().pullPlugins(['myorg/nf-foo@1.0.0'])
        assert store.toFile().list() as List == ['myorg~nf-foo-1.0.0']
        and:
        registry(entry('myorg/nf-foo', 'myorg/nf-foo', z))

        when:
        updater().pullPlugins(['myorg/nf-foo@1.0.0'])

        then:
        registryCalls() == 2
        isLink('nf-foo-1.0.0', 'myorg~nf-foo-1.0.0')
    }
}
