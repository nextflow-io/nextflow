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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.CompoundPluginRepository
import org.pf4j.PluginDescriptorFinder
import org.pf4j.PluginLoader
import org.pf4j.PluginRepository
/**
 * Plugin manager specialised for unit testing.
 *
 * The plugin root must be defined the sub-project root directory
 * e.g. 'ROOT/modules/nf-commons'.
 *
 * The plugin code must be defined in the 'testFixtures' source tree
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TestPluginManager extends DevPluginManager {

    private Path pluginRoot

    TestPluginManager(Path root) {
        super(root)
        this.pluginRoot = root
    }

    @Override
    protected PluginDescriptorFinder createPluginDescriptorFinder() {
        return new TestPluginDescriptorFinder()
    }

    @Override
    protected PluginLoader createPluginLoader() {
        return new TestPluginLoader(this)
    }

    @Override
    protected PluginRepository createPluginRepository() {
        def repos = new CompoundPluginRepository()
        // main dev repo
        final root = getPluginsRoot()
        log.debug "Added plugin root repository: ${root}"
        repos.add( new PluginRepository() {
            @Override
            List<Path> getPluginPaths() {
                return List.of(pluginRoot)
            }

            @Override
            boolean deletePluginPath(Path pluginPath) {
                log.debug "Test mode -- Ignore deletePluginPath('$pluginPath')"
                return false
            }
        })

        return repos
    }
}
