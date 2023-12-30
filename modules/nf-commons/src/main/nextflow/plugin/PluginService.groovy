/*
 * Copyright 2020-2022, Seqera Labs
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

import org.pf4j.PluginManager

/**
 * Define the contract for plugin service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface PluginService {

    void init()

    void init(boolean embeddedMode)

    void init(Path root, String mode, CustomPluginManager pluginManager)

    void load(Map config)

    void setup()

    void setup(Map config)

    void start(String pluginId)

    void stop()

    <T> List<T> getExtensions(Class<T> type)

    <T> List<T> getExtensions(Class<T> type, String pluginId)

    <T> List<T> getPriorityExtensions(Class<T> type)

    <T> List<T> getPriorityExtensions(Class<T> type, String group)

    void pullPlugins(List<String> ids)

    boolean startIfMissing(String pluginId)

    boolean isStarted(String pluginId)

    PluginManager getManager()

}
