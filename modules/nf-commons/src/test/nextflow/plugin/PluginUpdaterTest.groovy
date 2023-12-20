package nextflow.plugin

import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

import nextflow.BuildInfo
import org.pf4j.Plugin
import org.pf4j.PluginDescriptor
import org.pf4j.PluginWrapper
import org.pf4j.update.PluginInfo
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginUpdaterTest extends Specification {

    static class FooPlugin extends Plugin {
        FooPlugin(PluginWrapper wrapper) {
            super(wrapper)
        }
    }


    def 'should install a plugin' () {
        given:
        def PLUGIN = 'my-plugin-1.0.0'
        def folder = Files.createTempDirectory('test')

        and:
        // the plugin to be installed
        Path plugin = folder.resolve('plugins')
        def plugin1 = createPlugin(plugin, 'my-plugin', '1.0.0', FooPlugin.class)
        def plugin2 = createPlugin(plugin, 'my-plugin', '2.0.0', FooPlugin.class)
        def zip1 = zipDir(plugin1)
        def zip2 = zipDir(plugin2)
        and:
        // this represents the remote repo from where plugins are downloaded
        def repoDir = Files.createDirectory(folder.resolve('repo'))
        createRepositoryIndex(repoDir, zip1, zip2)

        and:
        // the central cache where downloaded unzipped plugins are kept
        def cacheDir = Files.createDirectory(folder.resolve('cache'))
        and:
        def manager = new LocalPluginManager(cacheDir)
        def updater = new PluginUpdater(manager, cacheDir, new URL("file:${repoDir.resolve('plugins.json')}"))

        when:
        updater.installPlugin( 'my-plugin', '1.0.0' )

        then:
        manager.getPlugin('my-plugin').plugin.class == FooPlugin.class
        manager.getPlugin('my-plugin').descriptor.getPluginId() == 'my-plugin'
        manager.getPlugin('my-plugin').descriptor.getVersion() == '1.0.0'
        and:
        cacheDir.resolve(PLUGIN).exists()
        cacheDir.resolve(PLUGIN).isDirectory()
        cacheDir.resolve(PLUGIN).resolve('MANIFEST.MF').isFile()
        and:
        manager.localRoot.resolve(PLUGIN).exists()
        manager.localRoot.resolve(PLUGIN).isLink()
        manager.localRoot.resolve(PLUGIN).resolve('MANIFEST.MF').text == cacheDir.resolve(PLUGIN).resolve('MANIFEST.MF').text

        cleanup:
        folder?.deleteDir()
    }


    def 'should update a plugin' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        // the plugin to be installed
        def pluginDir = folder.resolve('plugins')
        def plugin1 = createPlugin(pluginDir, 'my-plugin', '1.0.0', FooPlugin.class)
        def plugin2 = createPlugin(pluginDir, 'my-plugin', '2.0.0', FooPlugin.class)
        def zip1 = zipDir(plugin1)
        def zip2 = zipDir(plugin2)
        and:
        // this represents the remote repo from where plugins are downloaded
        def repoDir = Files.createDirectory(folder.resolve('repo'))
        createRepositoryIndex(repoDir, zip1, zip2)
        and:
        // the central cache where downloaded unzipped plugins are kept
        def cacheDir = Files.createDirectory(folder.resolve('cache'))
        and:
        def manager = new LocalPluginManager(cacheDir)
        def updater = new PluginUpdater(manager, cacheDir, new URL("file:${repoDir.resolve('plugins.json')}"))

        when:
        updater.installPlugin('my-plugin', '1.0.0')
        then:
        manager.getPlugin('my-plugin').plugin.class == FooPlugin.class
        manager.getPlugin('my-plugin').descriptor.getPluginId() == 'my-plugin'
        manager.getPlugin('my-plugin').descriptor.getVersion() == '1.0.0'
        and:
        cacheDir.resolve('my-plugin-1.0.0').exists()
        cacheDir.resolve('my-plugin-1.0.0').isDirectory()

        when:
        updater.updatePlugin( 'my-plugin', '2.0.0' )
        then:
        manager.localRoot.resolve('my-plugin-2.0.0').exists()
        !manager.localRoot.resolve('my-plugin-1.0.0').exists()
        and:
        cacheDir.resolve('my-plugin-1.0.0').exists()
        cacheDir.resolve('my-plugin-2.0.0').exists()
        and:
        manager.getPlugin('my-plugin').plugin.class == FooPlugin.class
        manager.getPlugin('my-plugin').descriptor.getPluginId() == 'my-plugin'
        manager.getPlugin('my-plugin').descriptor.getVersion() == '2.0.0'

        cleanup:
        folder?.deleteDir()
    }


    def 'should not download existing plugin' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        // this represents the remote repo from where plugins are downloaded
        def repoDir = Files.createDirectory(folder.resolve('repo'))
        createEmptyIndex(repoDir)
        and:
        // the central cache where downloaded unzipped plugins are kept
        def cacheDir = Files.createDirectory(folder.resolve('cache'))
        def plugin1 = createPlugin(cacheDir,'my-plugin', '1.0.0', FooPlugin.class)
        def plugin2 = createPlugin(cacheDir,'my-plugin', '2.0.0', FooPlugin.class)
        and:
        def manager = new LocalPluginManager(cacheDir)
        def updater = new PluginUpdater(manager, cacheDir, new URL("file:${repoDir.resolve('plugins.json')}"))

        when:
        updater.installPlugin('my-plugin', '1.0.0')
        then:
        manager.getPlugin('my-plugin').plugin.class == FooPlugin.class
        manager.getPlugin('my-plugin').descriptor.getPluginId() == 'my-plugin'
        manager.getPlugin('my-plugin').descriptor.getVersion() == '1.0.0'
        and:
        cacheDir.resolve('my-plugin-1.0.0').exists()
        cacheDir.resolve('my-plugin-1.0.0').isDirectory()
        and:
        manager.localRoot.resolve('my-plugin-1.0.0').exists()
        manager.localRoot.resolve('my-plugin-1.0.0').isLink()

        when:
        updater.updatePlugin( 'my-plugin', '2.0.0' )
        then:
        manager.localRoot.resolve('my-plugin-2.0.0').exists()
        !manager.localRoot.resolve('my-plugin-1.0.0').exists()
        and:
        cacheDir.resolve('my-plugin-1.0.0').exists()
        cacheDir.resolve('my-plugin-2.0.0').exists()
        and:
        manager.getPlugin('my-plugin').plugin.class == FooPlugin.class
        manager.getPlugin('my-plugin').descriptor.getPluginId() == 'my-plugin'
        manager.getPlugin('my-plugin').descriptor.getVersion() == '2.0.0'

        cleanup:
        folder?.deleteDir()
    }


    static private Path createPlugin(Path baseDir, String id, String ver, Class clazz) {
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
                """
                .stripIndent()

        return pluginDir
    }

    static private void createRepositoryIndex(Path repoDir, Path zip1, Path zip2) {
        repoDir.resolve('plugins.json').text = """
              [{
                "id": "my-plugin",
                "description": "Test plugin",
                "releases": [
                  {
                    "version": "1.0.0",
                    "date": "Jun 25, 2020 9:58:35 PM",
                    "url": "file:${zip1}"
                  },
                  {
                    "version": "2.0.0",
                    "date": "Jun 25, 2020 9:58:35 PM",
                    "url": "file:${zip2}"
                  }
                ]
              }]
            """
    }

    static private void createEmptyIndex(Path repoDir) {
        repoDir.resolve('plugins.json').text = """
              [{
                "id": "my-plugin",
                "description": "Test plugin",
                "releases": [ ]
              }]
            """
    }

    static Path zipDir(final Path folder) throws IOException {

        def zipFilePath = folder.resolveSibling( "${folder.name}.zip" )

        try (
                FileOutputStream fos = new FileOutputStream(zipFilePath.toFile());
                ZipOutputStream zos = new ZipOutputStream(fos)
        ) {
            Files.walkFileTree(folder, new SimpleFileVisitor<Path>() {
                public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
                    zos.putNextEntry(new ZipEntry(folder.relativize(file).toString()));
                    Files.copy(file, zos);
                    zos.closeEntry();
                    return FileVisitResult.CONTINUE;
                }

                public FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException {
                    zos.putNextEntry(new ZipEntry(folder.relativize(dir).toString() + "/"));
                    zos.closeEntry();
                    return FileVisitResult.CONTINUE;
                }
            });
        }

        return zipFilePath
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

    def 'should save move plugin directory' () {
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
}
