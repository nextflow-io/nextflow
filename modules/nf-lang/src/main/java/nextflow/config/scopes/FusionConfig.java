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
import nextflow.script.types.MemoryUnit;

public class FusionConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable/disable the use of Fusion file system.
    """)
    public boolean enabled;

    @ConfigOption
    @Description("""
        The maximum size of the local cache used by the Fusion client.
    """)
    public MemoryUnit cacheSize;

    @ConfigOption
    @Description("""
        The URL from where the container layer provisioning the Fusion client is downloaded.
    """)
    public String containerConfigUrl;

    @ConfigOption
    @Description("""
        When `true` the access credentials required by the underlying object storage are exported to the task execution environment.
    """)
    public boolean exportStorageCredentials;

    @ConfigOption
    @Description("""
        The level of logging emitted by the Fusion client.
    """)
    public String logLevel;

    @ConfigOption
    @Description("""
        Where the logging output is written. 
    """)
    public String logOutput;

    @ConfigOption
    @Description("""
        Enables the use of privileged containers when using Fusion (default: `true`).
    """)
    public boolean privileged;

    @ConfigOption
    @Description("""
        The pattern that determines how tags are applied to files created via the Fusion client (default: `[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)`). Set to `false` to disable tags.
    """)
    public String tags;

}
