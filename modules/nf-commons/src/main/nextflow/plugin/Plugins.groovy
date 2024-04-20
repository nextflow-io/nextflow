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
import org.pf4j.PluginManager
/**
 * Plugin manager specialized for Nextflow build environment
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Plugins {

    public static final String DEFAULT_PLUGINS_REPO = 'https://raw.githubusercontent.com/nextflow-io/plugins/main/plugins.json'

    private final static PluginsFacade INSTANCE = new PluginsFacade()

    static PluginManager getManager() { INSTANCE.manager }

    static synchronized void init(boolean embeddedMode=false) {
        INSTANCE.init(embeddedMode)
    }

    static synchronized void init(Path root, String mode, CustomPluginManager pluginManager) {
        INSTANCE.init(root, mode, pluginManager)
    }

    /**
     * @param config
     */
    static void load(Map config) {
        INSTANCE.load(config)
    }

    static void start(String pluginId) {
        INSTANCE.start(pluginId)
    }

    static synchronized void stop() {
        INSTANCE.stop()
    }

    static <T> List<T> getExtensions(Class<T> type) {
        INSTANCE.getExtensions(type)
    }

    static <T> List<T> getPriorityExtensions(Class<T> type, String group=null) {
        INSTANCE.getPriorityExtensions(type,group)
    }

    static <T> T getExtension(Class<T> type) {
        final allExtensions = INSTANCE.getExtensions(type)
        return allExtensions ? allExtensions.first() : null
    }

    static <T> List<T> getExtensionsInPluginId(Class<T> type, String pluginId) {
        final allExtensions = INSTANCE.getExtensions(type, pluginId)
        return allExtensions
    }

    static void pull(List<String> ids) {
        INSTANCE.pullPlugins(ids)
    }

    static boolean startIfMissing(String pluginId) {
        if( INSTANCE ) {
            return INSTANCE.startIfMissing(pluginId)
        } else {
            log.debug "Plugins subsystem not available - Ignoring installIfMissing('$pluginId')"
            return false
        }
    }

    static boolean isStarted(String pluginId) {
        INSTANCE.isStarted(pluginId)
    }
}
