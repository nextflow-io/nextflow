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

public class AwsConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        AWS account access key.
    """)
    public String accessKey;

    @ConfigOption
    @Description("""
        AWS profile from `~/.aws/credentials`.
    """)
    public String profile;

    @ConfigOption
    @Description("""
        AWS region (e.g. `us-east-1`).
    """)
    public String region;

    @ConfigOption
    @Description("""
        AWS account secret key.
    """)
    public String secretKey;

    public AwsBatchConfig batch;

    public AwsClientConfig client;

}
