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
        return new ManifestPluginDescriptorFinder()
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

    /**
     * Accept only the plugin directories that have been built, that is those holding
     * a manifest in one of the dev classpath directories.
     *
     * The manifest locations are derived from {@link DevPluginClasspath} so they stay
     * in sync with the directories used by the dev plugin class loader.
     */
    @CompileStatic
    private static class BuiltPluginFilter implements FileFilter {

        private static final Collection<String> CLASSES_DIRS = new DevPluginClasspath().getClassesDirectories()

        @Override
        boolean accept(File file) {
            if( !file.isDirectory() || file.isHidden() )
                return false
            return CLASSES_DIRS.any { new File(file, "$it/META-INF/MANIFEST.MF").isFile() }
        }
    }

}
