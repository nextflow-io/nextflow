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

package nextflow.util

import static nextflow.util.ConfigHelper.*

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
/**
 * Models retry policy configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class RetryConfig {

    private final static Duration DEFAULT_DELAY = Duration.of('350ms')
    private final static Duration DEFAULT_MAX_DELAY = Duration.of('90s')
    private final static Integer DEFAULT_MAX_ATTEMPTS = 5
    private final static Double DEFAULT_JITTER = 0.25

    private final static String ENV_PREFIX = 'NXF_RETRY_POLICY_'

    final private Duration delay
    final private Duration maxDelay
    final private int maxAttempts
    final private double jitter

    RetryConfig() {
        this(Collections.emptyMap())
    }

    RetryConfig(Map config) {
        delay =
            valueOf(config, 'delay', ENV_PREFIX, DEFAULT_DELAY, Duration)
        maxDelay =
            valueOf(config, 'maxDelay', ENV_PREFIX, DEFAULT_MAX_DELAY, Duration)
        maxAttempts =
            valueOf(config, 'maxAttempts', ENV_PREFIX, DEFAULT_MAX_ATTEMPTS, Integer)
        jitter =
            valueOf(config, 'jitter', ENV_PREFIX, DEFAULT_JITTER, Double)
    }

    Duration getDelay() { delay }

    Duration getMaxDelay() { maxDelay }

    int getMaxAttempts() { maxAttempts }

    double getJitter() { jitter }

    static RetryConfig config() {
        config(Global.session as Session)
    }

    static RetryConfig config(Session session) {
        if( session ) {
            return new RetryConfig(session.config.navigate('nextflow.retryPolicy') as Map ?: Collections.emptyMap())
        }
        log.warn "Missing nextflow session - using default retry config"
        return new RetryConfig()
    }


}
