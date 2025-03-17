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
package nextflow.config.scopes;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.script.dsl.Description;
import nextflow.script.types.Duration;

public class WaveRetryConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The initial delay when a failing HTTP request is retried (default: `150ms`).
    """)
    public Duration delay;

    @ConfigOption
    @Description("""
        The jitter factor used to randomly vary retry delays (default: `0.25`).
    """)
    public double jitter;

    @ConfigOption
    @Description("""
        The max number of attempts a failing HTTP request is retried (default: `5`).
    """)
    public int maxAttempts;

    @ConfigOption
    @Description("""
        The max delay when a failing HTTP request is retried (default: `90 seconds`).
    """)
    public Duration maxDelay;

}
