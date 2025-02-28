package nextflow.plugin

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

import nextflow.BuildInfo
import com.github.zafarkhaja.semver.Version
import org.pf4j.Plugin
import org.pf4j.PluginDescriptor
import org.pf4j.PluginRuntimeException
import org.pf4j.PluginWrapper
import org.pf4j.update.PluginInfo
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginUpdaterTest extends Specification {

    private static final String PLUGIN_ID = 'my-plugin'

    static class FooPlugin extends Plugin {
        FooPlugin(PluginWrapper wrapper) {
            super(wrapper)
        }
    }

    static class MockPlugin {
        String version
        Path path
        Path zip
    }

    // ---------------------------

    def 'should install a plugin' () {
        given:
        def PLUGIN = 'my-plugin-1.0.0'
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.0.0', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), [])
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, false)

        when:
        updater.installPlugin( PLUGIN_ID, '1.0.0' )

        then:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '1.0.0'
        and:
        local.resolve(PLUGIN).exists()
        local.resolve(PLUGIN).isDirectory()
        local.resolve(PLUGIN).resolve('MANIFEST.MF').isFile()
        and:
        manager.localRoot.resolve(PLUGIN).exists()
        manager.localRoot.resolve(PLUGIN).isLink()
        manager.localRoot.resolve(PLUGIN).resolve('MANIFEST.MF').text == local.resolve(PLUGIN).resolve('MANIFEST.MF').text

        cleanup:
        folder?.deleteDir()
    }


    def 'should update a plugin' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.0.0', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), [])
        and:
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, false)

        when:
        updater.installPlugin(PLUGIN_ID, '1.0.0')
        then:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '1.0.0'
        and:
        local.resolve('my-plugin-1.0.0').exists()
        local.resolve('my-plugin-1.0.0').isDirectory()

        when:
        updater.updatePlugin( PLUGIN_ID, '2.0.0' )
        then:
        manager.localRoot.resolve('my-plugin-2.0.0').exists()
        !manager.localRoot.resolve('my-plugin-1.0.0').exists()
        and:
        local.resolve('my-plugin-1.0.0').exists()
        local.resolve('my-plugin-2.0.0').exists()
        and:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '2.0.0'

        cleanup:
        folder?.deleteDir()
    }


    def 'should not download existing plugin' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.0.0', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), [])
        and:
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, false)

        when:
        updater.installPlugin(PLUGIN_ID, '1.0.0')
        then:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '1.0.0'
        and:
        local.resolve('my-plugin-1.0.0').exists()
        local.resolve('my-plugin-1.0.0').isDirectory()
        and:
        manager.localRoot.resolve('my-plugin-1.0.0').exists()
        manager.localRoot.resolve('my-plugin-1.0.0').isLink()

        when:
        updater.updatePlugin( PLUGIN_ID, '2.0.0' )
        then:
        manager.localRoot.resolve('my-plugin-2.0.0').exists()
        !manager.localRoot.resolve('my-plugin-1.0.0').exists()
        and:
        local.resolve('my-plugin-1.0.0').exists()
        local.resolve('my-plugin-2.0.0').exists()
        and:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '2.0.0'

        cleanup:
        folder?.deleteDir()
    }


    def 'resolve plugin version from range' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.2.0', '1.2.3', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), [])
        and:
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, false)

        when:
        def success = updater.installPlugin(PLUGIN_ID, '~1.2.0')
        then:
        success
        and:
        manager.getPlugin(PLUGIN_ID).plugin.class == FooPlugin.class
        manager.getPlugin(PLUGIN_ID).descriptor.getPluginId() == PLUGIN_ID
        manager.getPlugin(PLUGIN_ID).descriptor.getVersion() == '1.2.3'
        and:
        local.resolve('my-plugin-1.2.3').exists()
        local.resolve('my-plugin-1.2.3').isDirectory()
        local.resolve('my-plugin-1.2.3').resolve('MANIFEST.MF').isFile()

        cleanup:
        folder?.deleteDir()
    }


    def 'in offline mode, should use already downloaded plugin' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.2.0', '1.2.3', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), ['1.2.1', '1.2.2'])
        and:
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, true)

        when:
        // 1.2.1 exists in local plugin cache
        def success = updater.installPlugin(PLUGIN_ID, '1.2.1')
        then:
        success
        manager.getPlugin(PLUGIN_ID).descriptor.version == '1.2.1'

        when:
        // 2.0.0 only exists in remote repo, should shouldn't resolve
        manager.unloadPlugin(PLUGIN_ID)
        updater.installPlugin(PLUGIN_ID, '2.0.0')
        then:
        def error = thrown(PluginRuntimeException.class)
        error.message == 'Plugin my-plugin with version @2.0.0 does not exist in the repository'

        cleanup:
        folder?.deleteDir()
    }


    def 'in offline mode, resolve plugin version from range' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def remote = remoteRepository(folder.resolve('repo'), ['1.2.0', '1.2.3', '2.0.0'])
        and:
        def local = localCache(folder.resolve('plugins'), ['1.2.1', '1.2.2'])
        and:
        def manager = new LocalPluginManager(local)
        def updater = new PluginUpdater(manager, local, remote, true)

        when:
        // should be able to find a match in range >=1.2.1 && <2.0.0
        def success = updater.installPlugin(PLUGIN_ID, '~1.2.1')
        then:
        success
        // version should be 1.2.2 (local) not 1.2.3 (remote)
        manager.getPlugin(PLUGIN_ID).descriptor.version == '1.2.2'

        when:
        // should not be able to find a match in range >=2.0.0 && <3.0.0
        manager.unloadPlugin(PLUGIN_ID)
        updater.installPlugin(PLUGIN_ID, '~2.0.0')
        then:
        def error = thrown(IllegalStateException.class)
        error.message == 'Cannot find version of my-plugin plugin matching ~2.0.0'

        cleanup:
        folder?.deleteDir()
    }


    def 'versions with ~ should behave like +' () {
        // this test checks our understanding on ~ semver rules is correct
        when:
        def condition = "~1.2"
        then:
        Version.parse("1.2.0").satisfies(condition)
        Version.parse("1.2.1").satisfies(condition)
        Version.parse("1.2.2").satisfies(condition)
        !Version.parse("1.3.0").satisfies(condition)

        when:
        condition = "~1.2.1"
        then:
        !Version.parse("1.2.0").satisfies(condition)
        Version.parse("1.2.1").satisfies(condition)
        Version.parse("1.2.2").satisfies(condition)
        !Version.parse("1.3.0").satisfies(condition)
    }


    def 'should find matching plugin' () {
        given:
        def r1 = new PluginInfo.PluginRelease(version: '1.4.0', url: 'http://xyz')
        def r2 = new PluginInfo.PluginRelease(version: '1.5.0', url: 'http://xyz', requires: ">=${BuildInfo.version}")
        def r3 = new PluginInfo.PluginRelease(version: '1.5.1', url: 'http://xyz', requires: ">=${BuildInfo.version}")
        def r4 = new PluginInfo.PluginRelease(version: '2.0.1', url: 'http://xyz', requires: ">=${BuildInfo.version}")
        def r5 = new PluginInfo.PluginRelease(version: '3.0.0', url: 'http://xyz', requires: '99.01.0')
        def PLUGINS = [
                'nf-foo': new PluginInfo(id:'nf-foo', releases: [r1, r2, r3, r4, r5]),
                'nf-bar': new PluginInfo(id:'nf-bar', releases: [])
        ]
        def manager = Mock(CustomPluginManager)
        PluginUpdater updater = Spy(PluginUpdater, constructorArgs: [manager])


        when:
        def ret = updater.findNewestMatchingRelease('nf-foo', '1.5.0')
        then:
        manager.getVersionManager() >> new CustomVersionManager()
        updater.getPluginsMap() >> PLUGINS
        and:
        ret == r2

        when:
        ret = updater.findNewestMatchingRelease('nf-foo', '1.5.*')
        then:
        manager.getVersionManager() >> new CustomVersionManager()
        updater.getPluginsMap() >> PLUGINS
        and:
        ret == r3

        when:
        ret = updater.findNewestMatchingRelease('nf-foo', '>=2.0')
        then:
        manager.getVersionManager() >> new CustomVersionManager()
        updater.getPluginsMap() >> PLUGINS
        and:
        ret == r4
    }


    def 'should safely move plugin directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        // create a plugin dir
        def source = Files.createDirectory(folder.resolve('plugin-1.0.0'))
        source.resolve('manifest.txt').text = 'Plugin-1.0.0'
        source.resolve('foo.txt').text = "I'm foo"
        source.resolve('bar.txt').text = "I'm bar"
        Files.createDirectory(source.resolve('alpha'))
        Files.createDirectory(source.resolve('alpha/beta'))
        source.resolve('alpha/one.txt').text = 'File 1'
        source.resolve('alpha/beta/two.txt').text = 'File 2'
        source.resolve('alpha/beta/three.txt').text = 'File 3'
        and:
        def target = folder.resolve('new-plugin-1.0')
        and:
        Files.createDirectory(target)
        target.resolve('foo').text = 'some content'
        and:
        def updater = new PluginUpdater(Mock(CustomPluginManager))

        when:
        updater.safeMove0(source, target)

        then:
        Files.exists(target)
        and:
        target.resolve('manifest.txt').text == 'Plugin-1.0.0'
        target.resolve('alpha/one.txt').text == 'File 1'
        target.resolve('alpha/beta/three.txt').text == 'File 3'
        and:
        !Files.exists(source)

        cleanup:
        folder.deleteDir()
    }


    def 'should move plugin directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        // create a plugin dir
        def source = Files.createDirectory(folder.resolve('plugin-1.0.0'))
        source.resolve('manifest.txt').text = 'Plugin-1.0.0'
        source.resolve('foo.txt').text = "I'm foo"
        source.resolve('bar.txt').text = "I'm bar"
        Files.createDirectory(source.resolve('alpha'))
        Files.createDirectory(source.resolve('alpha/beta'))
        source.resolve('alpha/one.txt').text = 'File 1'
        source.resolve('alpha/beta/two.txt').text = 'File 2'
        source.resolve('alpha/beta/three.txt').text = 'File 3'
        and:
        def target = folder.resolve('new-plugin-1.0')
        and:
        Files.createDirectory(target); target.resolve('foo').text = 'some content'
        and:
        def updater = new PluginUpdater(Mock(CustomPluginManager))

        when:
        updater.safeMove(source, target)

        then:
        Files.exists(target)
        and:
        target.resolve('manifest.txt').text == 'Plugin-1.0.0'
        target.resolve('alpha/one.txt').text == 'File 1'
        target.resolve('alpha/beta/three.txt').text == 'File 3'
        and:
        !Files.exists(source)

        cleanup:
        folder.deleteDir()
    }


    @Unroll
    def 'validate should update method' () {
        given:
        def manager = Mock(CustomPluginManager) {
            getVersionManager() >> new CustomVersionManager()
        }
        and:
        def updater = new PluginUpdater(manager)

        and:
        def wrapper = Mock(PluginWrapper) {
            getDescriptor() >> Mock(PluginDescriptor) {
                getVersion() >> CURRENT
            }
        }

        expect:
        updater.shouldUpdate(ID, REQUIRED, wrapper) == EXPECT

        where:
        ID              | REQUIRED      | CURRENT | EXPECT
        'nf-amazon'     | '1.0.0'       | '1.0.0'   | false
        'nf-amazon'     | '1.0.0'       | '1.1.0'   | false
        'nf-amazon'     | '1.1.0'       | '1.0.0'   | true
    }


    @Unroll
    def 'should match meta file name' () {
        when:
        def matcher = PluginUpdater.META_REGEX.matcher(FILE_NAME)
        then:
        matcher.matches() == EXPECTED
        !EXPECTED || matcher.group(1) == PLUGIN
        
        where:
        FILE_NAME                               | EXPECTED  | PLUGIN
        'foo'                                   | false     | null
        'foo.json'                              | false     | null
        'nf-foo-1.0.0-meta.json'                | true      | 'nf-foo'
        'xpack-google-1.0.0-beta.3-meta.json'   | true      | 'xpack-google'
    }

    // -------------------------------------------------------------------------------------
    // setup helpers

    static private URL remoteRepository(Path dir, List<String> versions) {
        Files.createDirectory(dir)

        List<MockPlugin> plugins = versions
            .collect ( version -> createPlugin(dir, version) )
            .collect { plugin ->
                plugin.zip = zipDir(plugin.path)
                plugin
            }

        return createRepositoryIndex(dir, plugins).toUri().toURL()
    }

    static private Path localCache(Path dir, List<String> versions) {
        Files.createDirectory(dir)
        versions.each { version -> createPlugin(dir, version) }
        return dir
    }

    static private MockPlugin createPlugin(Path baseDir, String ver) {
        def id = "my-plugin"
        def clazz = FooPlugin.class
        def fqn = "$id-$ver".toString()
        def pluginDir = baseDir.resolve(fqn)
        pluginDir.mkdirs()

        pluginDir.resolve('file1.txt').text = 'foo'
        pluginDir.resolve('file2.txt').text = 'bar'
        pluginDir.resolve('MANIFEST.MF').text = """\
                Manifest-Version: 1.0
                Plugin-Class: ${clazz.getName()}
                Plugin-Id: $id
                Plugin-Version: $ver
                """.stripIndent()

        return new MockPlugin(version: ver, path: pluginDir)
    }

    static private Path createRepositoryIndex(Path repoDir, List<MockPlugin> plugins) {
        String releases = plugins
            .collect ( p ->
                """
                {
                "version": "${p.version}",
                "date": "Jun 25, 2020 9:58:35 PM",
                "url": "file:${p.zip}"
                }
                """
            )
            .join(",")

        Path index = repoDir.resolve('plugins.json')
        index.text = """
            [{
            "id": "my-plugin",
            "description": "Test plugin",
            "releases": [${releases}]
            }]
            """
        return index
    }

    static private Path zipDir(final Path folder) throws IOException {
        def zipFilePath = folder.resolveSibling( "${folder.name}.zip" )

        try (
            FileOutputStream fos = new FileOutputStream(zipFilePath.toFile());
            ZipOutputStream zos = new ZipOutputStream(fos)
        ) {
            Files.walkFileTree(folder, new SimpleFileVisitor<Path>() {
                FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    zos.putNextEntry(new ZipEntry(folder.relativize(file).toString()));
                    Files.copy(file, zos);
                    zos.closeEntry();
                    return FileVisitResult.CONTINUE;
                }

                FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                    zos.putNextEntry(new ZipEntry(folder.relativize(dir).toString() + "/"));
                    zos.closeEntry();
                    return FileVisitResult.CONTINUE;
                }
            });
        }

        return zipFilePath
    }
}
