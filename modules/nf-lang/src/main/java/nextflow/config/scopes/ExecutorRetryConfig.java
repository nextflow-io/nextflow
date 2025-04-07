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

public class ExecutorRetryConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Delay when retrying failed job submissions (default: `500ms`).
        """)
    public Duration delay;

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Jitter value when retrying failed job submissions (default: `0.25`).
        """)
    public double jitter;

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Max attempts when retrying failed job submissions (default: `3`).
        """)
    public int maxAttempts;

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Max delay when retrying failed job submissions (default: `30s`).
        """)
    public Duration maxDelay;

    @ConfigOption
    @Description("""
        *Used only by grid executors.*

        Regex pattern that when verified causes a failed submit operation to be re-tried (default: `Socket timed out`).
        """)
    public String reason;


}
