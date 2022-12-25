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

package io.seqera.wave.plugin.config

/**
 * Conda build options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondaOpts {

    final public String DEFAULT_MAMBA_IMAGE = 'mambaorg/micromamba:1.0.0'

    final String mambaImage
    final List<String> commands

    CondaOpts(Map opts) {
        this.mambaImage = opts.mambaImage ?: DEFAULT_MAMBA_IMAGE
        this.commands = opts.commands as List<String>
    }

}
