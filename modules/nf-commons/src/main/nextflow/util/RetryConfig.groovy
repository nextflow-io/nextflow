/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.util

import com.google.common.base.CaseFormat
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import io.seqera.util.retry.Retryable
import nextflow.Global
import nextflow.SysEnv
/**
 * Models retry policy configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class RetryConfig implements Retryable.Config {

    private final static Duration DEFAULT_DELAY = Duration.of('350ms')
    private final static Duration DEFAULT_MAX_DELAY = Duration.of('90s')
    private final static Integer DEFAULT_MAX_ATTEMPTS = 5
    private final static Double DEFAULT_JITTER = 0.25
    static final public double DEFAULT_MULTIPLIER = 2.0

    private final static String ENV_PREFIX = 'NXF_RETRY_POLICY_'

    final private Duration delay
    final private Duration maxDelay
    final private int maxAttempts
    final private double jitter
    final private double multiplier

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
        multiplier =
            valueOf(config, 'multiplier', ENV_PREFIX, DEFAULT_MULTIPLIER, Double)
    }

    Duration getDelay() {
        return delay
    }

    Duration getMaxDelay() {
        return maxDelay
    }

    @Override
    int getMaxAttempts() {
        return maxAttempts
    }

    @Override
    double getJitter() {
        return jitter
    }

    @Override
    double getMultiplier() {
        return multiplier
    }

    static RetryConfig config() {
        config(Global.config)
    }

    @Memoized
    static RetryConfig config(Map config) {
        if( config==null ) {
            log.debug "Missing nextflow session - using default retry config"
            config = Collections.emptyMap()
        }

        return new RetryConfig(getNestedConfig(config, 'nextflow', 'retryPolicy') ?: Collections.emptyMap())
    }

    private static Map getNestedConfig(Map config, String... keys) {
        def current = config
        for (String key : keys) {
            if (current instanceof Map && current.containsKey(key)) {
                current = current.get(key)
            } else {
                return null
            }
        }
        return current instanceof Map ? current : null
    }

    static <T> T valueOf(Map config, String name, String prefix, T defValue, Class<T> type)  {
        assert name, "Argument 'name' cannot be null or empty"
        assert type, "Argument 'type' cannot be null"

        // try to get the value from the config map
        final cfg = config?.get(name)
        if( cfg != null ) {
            return toType(cfg, type)
        }
        // try to fallback to the sys environment
        if( !prefix.endsWith('_') )
            prefix += '_'
        final key = prefix.toUpperCase() + CaseFormat.LOWER_CAMEL.to(CaseFormat.UPPER_UNDERSCORE, name)
        final env = SysEnv.get(key)
        if( env != null ) {
            return toType(env, type)
        }
        // return the default value
        return defValue
    }

    @CompileDynamic
    static protected <T> T toType(Object value, Class<T> type)  {
        if( value == null )
            return null
        if( type==Boolean.class ) {
            return type.cast(Boolean.valueOf(value.toString()))
        }
        else {
            return value.asType(type)
        }
    }
}
