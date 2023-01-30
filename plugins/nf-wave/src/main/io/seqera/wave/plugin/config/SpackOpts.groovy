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
 * Spack build options
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
class SpackOpts {

    final public String DEFAULT_SPACK_BUILDER_IMAGE = 'spack/ubuntu-jammy:v0.19.0'
    final public String DEFAULT_SPACK_RUNNER_IMAGE = 'ubuntu:22.04'

    final String spackBuilderImage
    final String spackRunnerImage
    final List<String> commands

    SpackOpts(Map opts) {
        this.spackBuilderImage = opts.spackBuilderImage ?: DEFAULT_SPACK_BUILDER_IMAGE
        this.spackRunnerImage = opts.spackRunnerImage ?: DEFAULT_SPACK_RUNNER_IMAGE
        this.commands = opts.commands as List<String>
    }

}
