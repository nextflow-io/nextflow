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
 */

package nextflow.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
/**
 * Represent Nextflow config as Map
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName('')
@CompileStatic
class ConfigMap extends LinkedHashMap implements ConfigScope {

    @ConfigOption
    @Description("""
        The remote work directory used by hybrid workflows. Equivalent to the `-bucket-dir` option of the `run` command.
    """)
    final String bucketDir

    @ConfigOption
    @Description("""
        Delete all files associated with a run in the work directory when the run completes successfully (default: `false`).
    """)
    final boolean cleanup

    @ConfigOption
    @Description("""
        The pipeline output directory. Equivalent to the `-output-dir` option of the `run` command.
    """)
    final String outputDir

    @ConfigOption
    @Description("""
        Enable the use of previously cached task executions. Equivalent to the `-resume` option of the `run` command.
    """)
    final boolean resume

    @ConfigOption
    @Description("""
        The pipeline work directory. Equivalent to the `-work-dir` option of the `run` command.
    """)
    final String workDir

    ConfigMap() {
    }

    ConfigMap(int initialCapacity) {
        super(initialCapacity)
    }

    ConfigMap(Map opts) {
        super(opts)
        bucketDir = opts.bucketDir
        cleanup = opts.cleanup as boolean
        outputDir = opts.outputDir
        resume = opts.resume as boolean
        workDir = opts.workDir
    }

}
