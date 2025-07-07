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

package nextflow.cloud.google.config

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 * Model Google storage retry settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GoogleRetryOpts {

    final int maxAttempts
    final double multiplier
    final Duration maxDelay

    GoogleRetryOpts(Map opts) {
        maxAttempts = opts.maxAttempts ? opts.maxAttempts as int : 10
        multiplier = opts.multiplier ? opts.multiplier as double : 2d
        maxDelay = opts.maxDelay ? opts.maxDelay as Duration : Duration.of('90s')
    }

    long maxDelaySecs() {
        return maxDelay.seconds
    }
}
