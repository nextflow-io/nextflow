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

import groovy.transform.CompileStatic

/**
 * Updater disabling plugin installation and update when
 * running the plugin manager in dev mode
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DevPluginUpdater extends PluginUpdater {

    DevPluginUpdater(CustomPluginManager pluginManager) {
        super(pluginManager)
    }

    @Override
    boolean installPlugin(String id, String version) {
        throw new UnsupportedOperationException("Install is not supported on dev mode - Missing plugin $id ${version?:''}")
    }

    @Override
    boolean updatePlugin(String id, String version) {
        throw new UnsupportedOperationException("Update is not supported on dev mode - Missing plugin $id ${version?:''}")
    }
}
