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

import java.nio.file.Path
import java.nio.file.Paths
import java.util.jar.Manifest

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.CompoundPluginRepository
import org.pf4j.DevelopmentPluginRepository
import org.pf4j.ManifestPluginDescriptorFinder
import org.pf4j.PluginDescriptorFinder
import org.pf4j.PluginLoader
import org.pf4j.PluginRepository
/**
 * Custom plugin manager that allow loading plugins from Groovy/Gradle/Intellij build environment
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class DevPluginManager extends CustomPluginManager {

    DevPluginManager(Path root) {
        super(root)
    }

    private List<Path> getExtensionRoots() {
        final ext = System .getenv('NXF_PLUGINS_DEV')
        if( !ext )
            return Collections.emptyList()

        return ext
                .tokenize(',')
                .collect { Paths.get(it).toAbsolutePath() }
    }

    @Override
    protected PluginDescriptorFinder createPluginDescriptorFinder() {
        return new DevManifestFinder()
    }

    @Override
    protected PluginLoader createPluginLoader() {
        return new DevPluginLoader(this)
    }

    @Override
    protected PluginRepository createPluginRepository() {
        def repos = new CompoundPluginRepository()
        // main dev repo
        log.debug "Add plugin root repository: ${getPluginsRoot()}"
        repos.add( devRepositoryFor(getPluginsRoot()) )
        // extension repos
        for( Path it : extensionRoots ) {
            log.debug "Add plugin dev repository: $it"
            repos.add( devRepositoryFor(it) )
        }
        return repos
    }

    private DevelopmentPluginRepository devRepositoryFor(Path root) {
        final repo = new DevelopmentPluginRepository(root)
        // Only consider plugin directories that have actually been compiled, i.e.
        // those for which a plugin manifest is available in the dev classpath. This
        // prevents pf4j from logging spurious "Cannot find the manifest path" errors
        // (with a full stack trace) for every plugin that has not been built yet in
        // the local dev environment - e.g. when running unit tests from a clean build.
        repo.setFilter(new BuiltPluginFilter())
        return repo
    }

    private static final Collection<String> CLASSES_DIRS = new DevPluginClasspath().getClassesDirectories()

    /**
     * Resolve the built plugin manifest that actually carries a {@code Plugin-Id}, looking only in the
     * dev classpath directories (derived from {@link DevPluginClasspath}) and in their fixed order.
     *
     * This is deterministic on purpose: a built plugin dir can hold more than one {@code MANIFEST.MF}
     * (e.g. a stray {@code build/tmp/jar/MANIFEST.MF} with no {@code Plugin-Id}), and pf4j's default
     * recursive lookup picks an arbitrary one depending on filesystem listing order - which
     * intermittently throws "Field 'id' cannot be empty". Returns null if none is found.
     */
    private static File builtManifest(File pluginDir) {
        for( String dir : CLASSES_DIRS ) {
            final f = new File(pluginDir, "$dir/META-INF/MANIFEST.MF")
            if( f.isFile() && f.withInputStream { new Manifest(it).mainAttributes.getValue('Plugin-Id') } )
                return f
        }
        return null
    }

    /**
     * Accept only the plugin directories that have been built, that is those holding
     * a valid (Plugin-Id bearing) manifest in one of the dev classpath directories.
     */
    @CompileStatic
    private static class BuiltPluginFilter implements FileFilter {
        @Override
        boolean accept(File file) {
            if( !file.isDirectory() || file.isHidden() )
                return false
            return builtManifest(file) != null
        }
    }

    /**
     * Read the plugin manifest from the deterministic dev-build location instead of pf4j's
     * default recursive first-match search over the whole plugin dir.
     */
    @CompileStatic
    private static class DevManifestFinder extends ManifestPluginDescriptorFinder {
        @Override
        protected Manifest readManifestFromDirectory(Path pluginPath) {
            final file = builtManifest(pluginPath.toFile())
            if( !file )
                return super.readManifestFromDirectory(pluginPath)
            return file.withInputStream { new Manifest(it) } as Manifest
        }
    }

}
