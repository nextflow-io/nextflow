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

package io.seqera.wave.plugin.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.trace.TraceHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
@Deprecated
class ReportOpts {

    final private Boolean enabled

    final private String file

    ReportOpts(Map opts) {
        this.enabled = opts.enabled as Boolean
        this.file = opts.file
    }

    boolean enabled() {
        enabled || file != null
    }

    String file() {
        if( file )
            return file
        return enabled
            ? "containers-${TraceHelper.launchTimestampFmt()}.config"
            : null
    }
}
