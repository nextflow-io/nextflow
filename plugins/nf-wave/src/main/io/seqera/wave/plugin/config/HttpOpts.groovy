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
import nextflow.util.Duration
import nextflow.util.RateUnit

/**
 * Model the HTTP client settings to connect the Wave service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class HttpOpts {

    final private Duration connectTimeout

    final private RateUnit maxRate

    HttpOpts(Map opts) {
        connectTimeout = opts.connectTimeout as Duration ?: Duration.of('30s')
        maxRate = opts.maxRate as RateUnit ?: RateUnit.of('1/sec')
    }

    java.time.Duration connectTimeout() {
        return java.time.Duration.ofMillis(connectTimeout.toMillis())
    }

    RateUnit maxRate() {
        return maxRate
    }
}
