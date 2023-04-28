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

    final public String DEFAULT_SPACK_BUILDER_IMAGE = 'spack/ubuntu-jammy:v0.19.2'
    final public String DEFAULT_SPACK_RUNNER_IMAGE = 'ubuntu:22.04'
    final public String DEFAULT_SPACK_OSPACKAGES = ''
    final public String DEFAULT_SPACK_FLAGS = '-O3'

    final Boolean checksum
    final String builderImage
    final String runnerImage
    final String osPackages
    final String cFlags
    final String cxxFlags
    final String fFlags
    final List<String> commands

    SpackOpts(Map opts) {
        this.checksum = ( opts.checksum != null ) ? opts.checksum : true
        this.builderImage = opts.builderImage ?: DEFAULT_SPACK_BUILDER_IMAGE
        this.runnerImage = opts.runnerImage ?: DEFAULT_SPACK_RUNNER_IMAGE
        this.osPackages = opts.osPackages ?: DEFAULT_SPACK_OSPACKAGES
        this.cFlags = opts.cFlags ?: DEFAULT_SPACK_FLAGS
        this.cxxFlags = opts.cxxFlags ?: DEFAULT_SPACK_FLAGS
        this.fFlags = opts.fFlags ?: DEFAULT_SPACK_FLAGS
        this.commands = opts.commands as List<String>
    }

}
