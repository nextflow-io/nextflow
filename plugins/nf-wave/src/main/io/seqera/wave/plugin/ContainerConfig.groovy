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

package io.seqera.wave.plugin

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.transform.builder.Builder

/**
 * Model a container configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Builder
@Canonical
@CompileStatic
@ToString(includePackage = false, includeNames = true)
class ContainerConfig {

    List<String> entrypoint
    List<String> cmd
    List<String> env
    String workingDir

    List<ContainerLayer> layers

    ContainerConfig appendLayer(ContainerLayer it)  {
        if( layers==null )
            layers = new ArrayList<>(10)
        layers.add(it)
        return this
    }

    ContainerConfig prependLayer(ContainerLayer it)  {
        if( layers==null )
            layers = new ArrayList<>(10)
        layers.add(0, it)
        return this
    }
}
