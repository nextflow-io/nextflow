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

public class AzureRegistryConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The container registry from which to pull the Docker images (default: `docker.io`).
    """)
    public String server;

    @ConfigOption
    @Description("""
        The username to connect to a private container registry.
    """)
    public String userName;

    @ConfigOption
    @Description("""
        The password to connect to a private container registry.
    """)
    public String password;

}
