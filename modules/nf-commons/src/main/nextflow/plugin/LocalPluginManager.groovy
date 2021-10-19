/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.plugin

import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.DefaultPluginLoader
import org.pf4j.DefaultPluginManager
import org.pf4j.ManifestPluginDescriptorFinder
import org.pf4j.PluginDescriptorFinder
import org.pf4j.PluginLoader
import org.pf4j.PluginRepository
import org.pf4j.PluginWrapper
/**
 * Custom plugin manager creating tracking plugins with symlinks to
 * the parent repository
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LocalPluginManager extends CustomPluginManager {

    private Path repository
    List<PluginSpec> specs

    LocalPluginManager(Path localRoot, Path repository, List<PluginSpec> specs) {
        super(localRoot)
        if( !localRoot ) throw new IllegalArgumentException("Missing Local plugins root directory")
        if( !repository ) throw new IllegalArgumentException("Missing plugins repository directory")
        this.repository = repository
        this.specs = specs
    }

    protected Path getLocalRoot() { getPluginsRoot() }

    @Override
    protected PluginDescriptorFinder createPluginDescriptorFinder() {
        return new ManifestPluginDescriptorFinder()
    }

    @Override
    protected PluginLoader createPluginLoader() {
        return new DefaultPluginLoader(this)
    }

    @Override
    protected PluginRepository createPluginRepository() {
        return new LocalPluginRepository(localRoot)
    }

    @Override
    void loadPlugins() {
        // sync local plugins
        for( PluginSpec it : specs ) { linkPlugin(it) }
        // proceed with loading
        super.loadPlugins()
    }

    protected void linkPlugin(PluginSpec spec) {
        if( spec.version ) {
            // given a plugin fully qualified, create a symlink from the repository to the local repo
            final name = "${spec.id}-${spec.version}"
            if( Files.exists(repository.resolve(name)) && !Files.exists(localRoot.resolve(name)) ) {
                final link = createLinkFromPath(repository.resolve(name))
                log.debug "Plugin $spec resolved to $link"
            }
        }
        else {
            final List<Path> avail = findPlugins(spec.id, repository)
            final List<Path> exist = findPlugins(spec.id, localRoot)
            if( avail && !exist ) {
                final link = createLinkFromPath(avail[0])
                log.debug "Plugin $spec resolved to $link"
            }
        }
    }

    protected List<Path> findPlugins(String name, Path repo) {
        final result = new ArrayList<Path>()
        repo.eachDir {
            if( it.fileName.toString().startsWith(name) ) {
                result.add(it)
            }
        }
        return result.sort()
    }

    /**
     * Override parent {@link DefaultPluginManager#loadPlugin(java.nio.file.Path)} method
     * creating a symlink the plugin path directory in the current plugin root
     *
     * @param pluginPath
     *      The file path where the plugin is stored unzipped) in the
     *      central plugins directory
     */
    @Override
    String loadPlugin(Path pluginPath) {
        final symlink = createLinkFromPath(pluginPath)
        super.loadPlugin(symlink)
    }

    @Override
    PluginWrapper loadPluginFromPath(Path pluginPath) {
        final symlink = createLinkFromPath(pluginPath)
        return super.loadPluginFromPath(symlink)
    }

    private Path createLinkFromPath(Path pluginPath) {
        if( pluginPath.startsWith(localRoot))
            return pluginPath
        if( !pluginPath )
            throw new IllegalArgumentException("Plugin path cannot be null")
        if( !Files.isDirectory(pluginPath) )
            throw new IllegalArgumentException("Plugin path must be a directory: $pluginPath")

        // create a symlink relative to the current root
        final symlink = localRoot.resolve(pluginPath.getFileName())
        createLink0(symlink, pluginPath)
        return symlink
    }

    private void createLink0(Path symlink, Path pluginPath) {
        try {
            log.trace "Creating local plugins root link: $symlink → $pluginPath"
            Files.createSymbolicLink(symlink, pluginPath)
        }
        catch (FileAlreadyExistsException e) {
            log.debug "Deleting existing local plugins root link: $symlink"
            if( !sameTarget(symlink, pluginPath) ) {
                Files.delete(symlink)
                Files.createSymbolicLink(symlink, pluginPath)
            }
        }
    }

    private boolean sameTarget(Path link, Path target) {
        try {
            return Files.readSymbolicLink(link) == target
        }
        catch (IOException e) {
            return false
        }
    }
}
