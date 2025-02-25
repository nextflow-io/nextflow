/*
 * Copyright 2013-2024, Seqera Labs
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
        repos.add( new DevelopmentPluginRepository(getPluginsRoot()) )
        // extension repos
        for( Path it : extensionRoots ) {
            log.debug "Add plugin dev repository: $it"
            repos.add( new DevelopmentPluginRepository(it) )
        }
        return repos
    }

}
