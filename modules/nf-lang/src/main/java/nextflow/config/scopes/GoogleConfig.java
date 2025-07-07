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

public class GoogleConfig implements ConfigScope {

    public GoogleBatchConfig batch;

    @ConfigOption
    @Description("""
        When `true`, the given Google Cloud project ID is used as the billing project for storage access (default: `false`). Required when accessing data from *requester pays enabled* buckets.

        [Read more](https://cloud.google.com/storage/docs/requester-pays)
    """)
    public boolean enableRequesterPaysBuckets;

    @ConfigOption
    @Description("""
        The HTTP connection timeout for Cloud Storage API requests (default: `'60s'`).
    """)
    public Duration httpConnectTimeout;

    @ConfigOption
    @Description("""
        The HTTP read timeout for Cloud Storage API requests (default: `'60s'`).
    """)
    public Duration httpReadTimeout;

    @ConfigOption
    @Description("""
        The Google Cloud location where jobs are executed (default: `us-central1`).
    """)
    public String location;

    @ConfigOption
    @Description("""
        The Google Cloud project ID to use for pipeline execution.
    """)
    public String project;

}
