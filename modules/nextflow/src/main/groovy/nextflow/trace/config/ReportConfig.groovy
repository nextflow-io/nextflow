/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.trace.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.trace.TraceHelper

@ScopeName("report")
@Description("""
    The `report` scope allows you to configure the workflow [execution report](https://nextflow.io/docs/latest/tracing.html#execution-report).
""")
@CompileStatic
class ReportConfig implements ConfigScope {

    static final int DEF_MAX_TASKS = 10_000

    @ConfigOption
    @Description("""
        Create the execution report on workflow completion (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The path of the created execution report file (default: `'report-<timestamp>.html'`).
    """)
    final String file

    @ConfigOption
    @Description("""
    """)
    final int maxTasks

    @ConfigOption
    @Description("""
        Overwrite any existing report file with the same name (default: `false`).
    """)
    final boolean overwrite

    /* required by extension point -- do not remove */
    ReportConfig() {}

    ReportConfig(Map opts) {
        enabled = opts.enabled as boolean
        file = opts.file ?: defaultFileName()
        maxTasks = opts.maxTasks != null ? opts.maxTasks as int : DEF_MAX_TASKS
        overwrite = opts.overwrite as boolean
    }

    static final String defaultFileName() {
        return "report-${TraceHelper.launchTimestampFmt()}.html"
    }

}
