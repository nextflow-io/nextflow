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

package io.seqera.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Top-level configuration scope for Seqera settings.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("seqera")
@Description("""
    The `seqera` scope provides configuration for Seqera services.
""")
@CompileStatic
class SeqeraConfig implements ConfigScope {

    @Description("""
        Configuration for the Seqera compute executor.
    """)
    final ExecutorOpts executor

    /* required by config scope -- do not remove */
    SeqeraConfig() {}

    SeqeraConfig(Map opts) {
        this.executor = opts.executor
            ? new ExecutorOpts(opts.executor as Map)
            : null
    }

    ExecutorOpts getExecutor() {
        return executor
    }
}
