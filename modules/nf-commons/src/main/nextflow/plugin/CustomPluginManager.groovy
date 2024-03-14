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
import nextflow.BuildInfo
import org.pf4j.DefaultPluginManager
import org.pf4j.ExtensionFactory
import org.pf4j.PluginWrapper
import org.pf4j.SingletonExtensionFactory
import org.pf4j.VersionManager
/**
 * Custom plugin manager to that allow accessing to {@code loadPluginFromPath} and
 * {@code resolvePlugins} method
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CustomPluginManager extends DefaultPluginManager {

    CustomPluginManager() {}

    CustomPluginManager(Path root) {
        super(root)
    }

    @Override
    PluginWrapper loadPluginFromPath(Path pluginPath) {
        super.loadPluginFromPath(pluginPath)
    }

    @Override
    void resolvePlugins() {
        super.resolvePlugins()
    }

    @Override
    protected VersionManager createVersionManager() {
        return new CustomVersionManager()
    }

    @Override
    String getSystemVersion() {
        return BuildInfo.version
    }

    @Override
    protected ExtensionFactory createExtensionFactory() {
        return new SingletonExtensionFactory(this)
    }
}
