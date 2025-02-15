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

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.pf4j.DefaultPluginRepository

/**
 * Extends the default plugin repository to avoid the deletion
 * of the plugin directory and delete instead the local symbolic link
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LocalPluginRepository extends DefaultPluginRepository {

    LocalPluginRepository(Path pluginsRoot) {
        super(pluginsRoot)
    }

    /**
     * Override {@link DefaultPluginRepository#deletePluginPath(java.nio.file.Path)}
     * to prevent deleting the real plugin directory
     *
     * @param pluginPath
     * @return
     */
    @Override
    boolean deletePluginPath(Path pluginPath) {
        if(Files.isSymbolicLink(pluginPath))
            try {
                Files.delete(pluginPath)
                return true
            }
            catch (Exception e) {
                log.debug "Unable to delete plugin path: $pluginPath"
                return false
            }
        else
            return super.deletePluginPath(pluginPath)
    }
}
