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

package nextflow.spack

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.Duration

/**
 * Model Spack configuration
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@ScopeName("spack")
@Description("""
    The `spack` scope controls the creation of a Spack environment by the Spack package manager.
""")
@CompileStatic
class SpackConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Execute tasks with Spack environments (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path where Spack environments are stored. It should be accessible from all compute nodes when using a shared file system.
    """)
    final String cacheDir

    @ConfigOption
    @Description("""
        Enable checksum verification of source tarballs (default: `true`).
    """)
    final boolean checksum

    @ConfigOption
    @Description("""
        The amount of time to wait for the Spack environment to be created before failing (default: `60 min`).
    """)
    final Duration createTimeout

    @ConfigOption
    @Description("""
        The maximum number of parallel package builds (default: the number of available CPUs).
    """)
    final Integer parallelBuilds

    /* required by extension point -- do not remove */
    SpackConfig() {}

    SpackConfig(Map opts, Map<String, String> env) {
        enabled = opts.enabled != null
            ? opts.enabled as boolean
            : (env.NXF_SPACK_ENABLED?.toString() == 'true')
        cacheDir = opts.cacheDir
        checksum = opts.checksum as boolean
        createTimeout = opts.createTimeout as Duration
        parallelBuilds = opts.parallelBuilds as Integer
    }
}
