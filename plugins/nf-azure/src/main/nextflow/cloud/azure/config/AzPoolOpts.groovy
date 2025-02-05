/*
 * Copyright 2021, Microsoft Corp
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
class AzPoolOpts implements CacheFunnel {

    static public final String DEFAULT_PUBLISHER = "microsoft-azure-batch"
    static public final String DEFAULT_OFFER = "ubuntu-server-container"
    static public final String DEFAULT_SKU = "batch.node.ubuntu 20.04"
    static public final String DEFAULT_VM_TYPE = "Standard_D4_v3"
    static public final OSType DEFAULT_OS_TYPE = OSType.LINUX
    static public final String DEFAULT_SHARE_ROOT_PATH = "/mnt/batch/tasks/fsmounts"
    static public final Duration DEFAULT_SCALE_INTERVAL = Duration.of('5 min')

    String runAs
    boolean privileged
    String publisher
    String offer
    String fileShareRootPath
    String sku
    OSType osType = DEFAULT_OS_TYPE
    ImageVerificationType verification = ImageVerificationType.VERIFIED

    String vmType
    Integer vmCount = 1
    boolean autoScale
    String scaleFormula
    Duration scaleInterval
    Integer maxVmCount

    String schedulePolicy // spread | pack
    String registry
    String userName
    String password

    String virtualNetwork
    boolean lowPriority
    AzStartTaskOpts startTask
    
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
