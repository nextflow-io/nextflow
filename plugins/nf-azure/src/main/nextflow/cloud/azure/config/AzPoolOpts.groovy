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

package nextflow.cloud.azure.config

import com.azure.compute.batch.models.ImageVerificationType
import com.azure.compute.batch.models.OSType
import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper
import nextflow.util.Duration
/**
 * Model the settings of a VM pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AzPoolOpts implements CacheFunnel, ConfigScope {

    static public final String DEFAULT_PUBLISHER = "microsoft-dsvm"
    static public final String DEFAULT_OFFER = "ubuntu-hpc"
    static public final String DEFAULT_SKU = "batch.node.ubuntu 24.04"
    static public final String DEFAULT_VM_TYPE = "Standard_D4a_v4"
    static public final OSType DEFAULT_OS_TYPE = OSType.LINUX
    static public final String DEFAULT_SHARE_ROOT_PATH = "/mnt/batch/tasks/fsmounts"
    static public final Duration DEFAULT_SCALE_INTERVAL = Duration.of('5 min')

    @ConfigOption
    @Description("""
        Enable autoscaling feature for the pool identified with `<name>`.
    """)
    final boolean autoScale

    @ConfigOption
    @Description("""
        The internal root mount point when mounting File Shares. Must be `/mnt/resource/batch/tasks/fsmounts` for CentOS nodes or `/mnt/batch/tasks/fsmounts` for Ubuntu nodes (default: CentOS).
    """)
    final String fileShareRootPath

    @ConfigOption
    @Description("""
        Enable the use of low-priority VMs (default: `false`).
    """)
    final boolean lowPriority

    @ConfigOption
    @Description("""
        The max number of virtual machines when using auto scaling.
    """)
    final Integer maxVmCount

    @ConfigOption
    @Description("""
        The mount options for mounting the file shares (default: `-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp`).
    """)
    final String mountOptions

    @ConfigOption
    @Description("""
        The offer type of the virtual machine type used by the pool identified with `<name>` (default: `centos-container`).
    """)
    final String offer

    @ConfigOption
    @Description("""
        Enable the task to run with elevated access. Ignored if `runAs` is set (default: `false`).
    """)
    final boolean privileged

    @ConfigOption
    @Description("""
        The publisher of virtual machine type used by the pool identified with `<name>` (default: `microsoft-azure-batch`).
    """)
    final String publisher

    @ConfigOption
    @Description("""
        The username under which the task is run. The user must already exist on each node of the pool.
    """)
    final String runAs

    @ConfigOption
    @Description("""
        The [scale formula](https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling) for the pool identified with `<name>`.
    """)
    final String scaleFormula

    @ConfigOption
    @Description("""
        The interval at which to automatically adjust the Pool size according to the autoscale formula. Must be at least 5 minutes and at most 168 hours (default: `10 mins`).
    """)
    final Duration scaleInterval

    @ConfigOption
    @Description("""
        The scheduling policy for the pool identified with `<name>`. Can be either `spread` or `pack` (default: `spread`).
    """)
    final String schedulePolicy

    @ConfigOption
    @Description("""
        The ID of the Compute Node agent SKU which the pool identified with `<name>` supports (default: `batch.node.centos 8`).
    """)
    final String sku

    final AzStartTaskOpts startTask

    @ConfigOption
    @Description("""
        The subnet ID of a virtual network in which to create the pool.
    """)
    final String virtualNetwork

    @ConfigOption
    @Description("""
        The number of virtual machines provisioned by the pool identified with `<name>`.
    """)
    final Integer vmCount

    @ConfigOption
    @Description("""
        The virtual machine type used by the pool identified with `<name>`.
    """)
    final String vmType

    OSType osType = DEFAULT_OS_TYPE
    ImageVerificationType verification = ImageVerificationType.VERIFIED

    String registry
    String userName
    String password

    AzPoolOpts() {
        this(Collections.emptyMap())
    }

    AzPoolOpts(Map opts) {
        this.runAs = opts.runAs ?: ''
        this.privileged = opts.privileged ?: false
        this.publisher = opts.publisher ?: DEFAULT_PUBLISHER
        this.offer = opts.offer ?: DEFAULT_OFFER
        this.sku = opts.sku ?: DEFAULT_SKU
        this.vmType = opts.vmType ?: DEFAULT_VM_TYPE
        this.fileShareRootPath = opts.fileShareRootPath ?: buildFileShareRootPath()
        this.vmCount = opts.vmCount as Integer ?: 1
        this.autoScale = opts.autoScale as boolean
        this.scaleFormula = opts.scaleFormula
        this.schedulePolicy = opts.schedulePolicy
        this.scaleInterval = opts.scaleInterval as Duration ?: DEFAULT_SCALE_INTERVAL
        this.maxVmCount = opts.maxVmCount as Integer ?: vmCount *3
        this.startTask = new AzStartTaskOpts( opts.startTask ? opts.startTask as Map : Map.of() )
        this.registry = opts.registry
        this.userName = opts.userName
        this.password = opts.password
        this.virtualNetwork = opts.virtualNetwork
        this.lowPriority = opts.lowPriority as boolean
    }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(runAs)
        hasher.putBoolean(privileged)
        hasher.putUnencodedChars(publisher)
        hasher.putUnencodedChars(offer)
        hasher.putUnencodedChars(sku)
        hasher.putUnencodedChars(vmType)
        hasher.putUnencodedChars(fileShareRootPath)
        hasher.putUnencodedChars(registry ?: '')
        hasher.putUnencodedChars(userName ?: '')
        hasher.putUnencodedChars(password ?: '')
        hasher.putInt(vmCount)
        hasher.putBoolean(autoScale)
        hasher.putUnencodedChars(scaleFormula ?: '')
        hasher.putUnencodedChars(schedulePolicy ?: '')
        hasher.putUnencodedChars(virtualNetwork ?: '')
        hasher.putBoolean(lowPriority)
        hasher.putUnencodedChars(startTask.script ?: '')
        hasher.putBoolean(startTask.privileged)
        return hasher
    }

    String buildFileShareRootPath() {
        if (this.sku ==~ /.*centos.*/)
            return "/mnt/resource/batch/tasks/fsmounts"
        else if (this.sku ==~ /.*ubuntu.*/)
            return "/mnt/batch/tasks/fsmounts"
        else
            return ''
	}
}
