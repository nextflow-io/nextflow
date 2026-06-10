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

package nextflow.container

import groovy.transform.CompileStatic

/**
 * Models generic container configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface ContainerConfig {

    default boolean canRunOciImage() {
        return false
    }

    default boolean entrypointOverride() {
        return ContainerHelper.entrypointOverride()
    }

    default String getFusionOptions() {
        return null
    }

    default Object getKill() {
        return null
    }

    default String getRegistry() {
        return null
    }

    default boolean getRegistryOverride() {
        return false
    }

    default String getTemp() {
        return null
    }

    String getEngine()

    List<String> getEnvWhitelist()

    boolean isEnabled()

}
