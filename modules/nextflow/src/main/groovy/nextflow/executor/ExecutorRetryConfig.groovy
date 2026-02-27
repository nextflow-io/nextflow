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
package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.Duration

@CompileStatic
class ExecutorRetryConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Delay when retrying failed job submissions (default: `500ms`).
    """)
    Duration delay = Duration.of('500ms')

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Jitter value when retrying failed job submissions (default: `0.25`).
    """)
    double jitter = 0.25

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Max attempts when retrying failed job submissions (default: `3`).
    """)
    int maxAttempts = 3

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Max delay when retrying failed job submissions (default: `30s`).
    """)
    Duration maxDelay = Duration.of('30s')

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Regex pattern that when verified causes a failed submit operation to be re-tried (default: `Socket timed out`).
    """)
    String reason = 'Socket timed out'

    ExecutorRetryConfig(Map opts) {
        if( opts.delay )
            delay = opts.delay as Duration
        if( opts.jitter )
            jitter = opts.jitter as double
        if( opts.maxAttempts )
            maxAttempts = opts.maxAttempts as int
        if( opts.maxDelay )
            maxDelay = opts.maxDelay as Duration
        if( opts.reason )
            reason = opts.reason
    }

}
