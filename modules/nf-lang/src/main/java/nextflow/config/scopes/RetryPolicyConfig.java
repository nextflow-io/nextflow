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
package nextflow.config.scopes;

import nextflow.config.spec.ConfigOption;
import nextflow.config.spec.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

/**
 * Models the `nextflow.retryPolicy` config scope.
 */
public class RetryPolicyConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Delay used for retryable operations (default: `350ms`).
    """)
    public Duration delay;

    @ConfigOption
    @Description("""
        Jitter value used for retryable operations (default: `0.25`).
    """)
    public Double jitter;

    @ConfigOption
    @Description("""
        Max attempts used for retryable operations (default: `5`).
    """)
    public Integer maxAttempts;

    @ConfigOption
    @Description("""
        Max delay used for retryable operations (default: `90s`).
    """)
    public Duration maxDelay;

    @ConfigOption
    @Description("""
        Multiplier used for the exponential backoff delay of retryable operations (default: `2.0`).
    """)
    public Double multiplier;

}
