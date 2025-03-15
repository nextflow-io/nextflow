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

import java.util.Map;

import nextflow.config.schema.ConfigOption;
import nextflow.config.schema.ConfigScope;
import nextflow.config.schema.PlaceholderName;
import nextflow.script.dsl.Description;

public class AzureBatchConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The batch service account name.
    """)
    public String accountName;

    @ConfigOption
    @Description("""
        The batch service account key.
    """)
    public String accountKey;

    @ConfigOption
    @Description("""
        Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: `false`).
    """)
    public boolean allowPoolCreation;

    @ConfigOption
    @Description("""
        Enable the automatic creation of batch pools depending on the pipeline resources demand (default: `true`).
    """)
    public String autoPoolMode;

    @ConfigOption
    @Description("""
        The mode in which the `azcopy` tool is installed by Nextflow (default: `'node'`). The following options are available:
        
        - `'node'`: the `azcopy` tool is installed once during the pool creation
        - `'task'`: the `azcopy` tool is installed for each task execution
        - `'off'`: the `azcopy` tool is not installed
    """)
    public String copyToolInstallMode;

    @ConfigOption
    @Description("""
        Delete all jobs when the workflow completes (default: `false`).
    """)
    public boolean deleteJobsOnCompletion;

    @ConfigOption
    @Description("""
        Delete all compute node pools when the workflow completes (default: `false`).
    """)
    public boolean deletePoolsOnCompletion;

    @ConfigOption
    @Description("""
        Delete each task when it completes (default: `true`).
    """)
    public boolean deleteTasksOnCompletion;

    @ConfigOption
    @Description("""
        The batch service endpoint e.g. `https://nfbatch1.westeurope.batch.azure.com`.
    """)
    public String endpoint;

    @ConfigOption
    @Description("""
        The name of the batch service region, e.g. `westeurope` or `eastus2`. Not needed when the endpoint is specified.
    """)
    public String location;

    @PlaceholderName("<name>")
    public Map<String, AzureBatchPoolConfig> pools;

    @ConfigOption
    @Description("""
        When the workflow completes, set all jobs to terminate on task completion (default: `true`).
    """)
    public boolean terminateJobsOnCompletion;

}
