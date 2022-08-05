package nextflow.plugin

import nextflow.NextflowMeta

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
trait BuildPluginTrait {

    Path buildPlugin(String name, String version, Path path) {
        def folder = Files.createTempDirectory('test')
        Path plugin = folder.resolve('plugins')
        def pluginPath = createPlugin(plugin, name, version, path)
        def zipPlugin = PluginUpdaterTest.zipDir(pluginPath)
        def repoDir = Files.createDirectory(folder.resolve('repo'))
        createRepositoryIndex(repoDir, zipPlugin)
        def cacheDir = Files.createDirectory(folder.resolve('cache'))
        def manager = new LocalPluginManager(cacheDir)
        Plugins.INSTANCE.manager = manager
        def updater = new PluginUpdater(manager, cacheDir, new URL("file:${repoDir.resolve('plugins.json')}"))
        updater.installPlugin( name, version )

        NextflowMeta.instance.strictMode(true)

        folder
    }


    private Path createPlugin(Path baseDir, String id, String ver, Path buildPath) {
        def fqn = "$id-$ver".toString()
        def pluginDir = baseDir.resolve(fqn)
        pluginDir.mkdirs()

        copyRecursive( Paths.get(buildPath.toString(), "classes/groovy/main"), pluginDir, Path.of('classes'))
        copyRecursive( Paths.get(buildPath.toString(), "classes/main"), pluginDir, Path.of('classes'))
        copyRecursive( Paths.get(buildPath.toString(), "classes/main"), pluginDir, Path.of('.'))

        pluginDir
    }

    private void copyRecursive(Path sourceDir, Path targetDir, Path subdir){
        // Traverse the file tree and copy each file/directory.
        Path finalDir = Paths.get(targetDir.toString(), subdir.toString())
        Files.walk(sourceDir)
                .forEach(sourcePath -> {
                    try {
                        Path targetPath = finalDir.resolve(sourceDir.relativize(sourcePath));
                        Files.copy(sourcePath, targetPath, StandardCopyOption.REPLACE_EXISTING);
                    } catch (IOException ex) {
                    }
                });
    }

    private void createRepositoryIndex(Path repoDir, Path zip) {
        repoDir.resolve('plugins.json').text = """
              [{
                "id": "nf-plugin-template",
                "description": "Test plugin",
                "releases": [
                  {
                    "version": "0.0.0",
                    "date": "Jun 25, 2020 9:58:35 PM",
                    "url": "file:${zip}"
                  }
                ]
              }]
            """
    }
}