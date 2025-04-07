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

public class AzureBatchPoolConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Enable autoscaling feature for the pool identified with `<name>`.
    """)
    public boolean autoScale;

    @ConfigOption
    @Description("""
        The internal root mount point when mounting File Shares. Must be `/mnt/resource/batch/tasks/fsmounts` for CentOS nodes or `/mnt/batch/tasks/fsmounts` for Ubuntu nodes (default: CentOS).
    """)
    public String fileShareRootPath;

    @ConfigOption
    @Description("""
        Enable the use of low-priority VMs (default: `false`).
    """)
    public boolean lowPriority;

    @ConfigOption
    @Description("""
        The max number of virtual machines when using auto scaling.
    """)
    public int maxVmCount;

    @ConfigOption
    @Description("""
        The mount options for mounting the file shares (default: `-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp`).
    """)
    public String mountOptions;

    @ConfigOption
    @Description("""
        The offer type of the virtual machine type used by the pool identified with `<name>` (default: `centos-container`).
    """)
    public String offer;

    @ConfigOption
    @Description("""
        Enable the task to run with elevated access. Ignored if `runAs` is set (default: `false`).
    """)
    public boolean privileged;

    @ConfigOption
    @Description("""
        The publisher of virtual machine type used by the pool identified with `<name>` (default: `microsoft-azure-batch`).
    """)
    public String publisher;

    @ConfigOption
    @Description("""
        The username under which the task is run. The user must already exist on each node of the pool.
    """)
    public String runAs;

    @ConfigOption
    @Description("""
        The scale formula for the pool identified with `<name>`.

        [Read more](https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling)
    """)
    public String scaleFormula;

    @ConfigOption
    @Description("""
        The interval at which to automatically adjust the Pool size according to the autoscale formula. Must be at least 5 minutes and at most 168 hours (default: `10 mins`).
    """)
    public Duration scaleInterval;

    @ConfigOption
    @Description("""
        The scheduling policy for the pool identified with `<name>`. Can be either `spread` or `pack` (default: `spread`).
    """)
    public String schedulePolicy;

    @ConfigOption
    @Description("""
        The ID of the Compute Node agent SKU which the pool identified with `<name>` supports (default: `batch.node.centos 8`).
    """)
    public String sku;

    public AzureBatchPoolStartTaskConfig startTask;

    @ConfigOption
    @Description("""
        The subnet ID of a virtual network in which to create the pool.
    """)
    public String virtualNetwork;

    @ConfigOption
    @Description("""
        The number of virtual machines provisioned by the pool identified with `<name>`.
    """)
    public String vmCount;

    @ConfigOption
    @Description("""
        The virtual machine type used by the pool identified with `<name>`.
    """)
    public String vmType;

}

class AzureBatchPoolStartTaskConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The `startTask` that is executed as the node joins the Azure Batch node pool.
    """)
    public String script;

    @ConfigOption
    @Description("""
        Enable the `startTask` to run with elevated access (default: `false`).
    """)
    public String privileged;

}
