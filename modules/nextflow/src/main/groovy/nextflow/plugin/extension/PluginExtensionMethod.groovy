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

package nextflow.plugin.extension

import groovy.transform.MapConstructor

/**
 * Hold a reference to a extension method provided by a Nextflow plugin
 */
@MapConstructor
class PluginExtensionMethod {
    /**
     * The name of the method that needs to be invoked
     */
    String method

    /**
     * The target object on which the method is going to be invoked
     */
    Object target

    /**
     * The plugin holding this extension
     */
    String pluginId
}
