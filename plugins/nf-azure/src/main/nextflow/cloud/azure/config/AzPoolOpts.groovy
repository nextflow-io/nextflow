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

import com.google.common.hash.Hasher
import com.microsoft.azure.batch.protocol.models.OSType
import com.microsoft.azure.batch.protocol.models.VerificationType
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper

/**
 * Model the setting of a VM pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AzPoolOpts implements CacheFunnel {

    static public final String DEFAULT_PUBLISHER = "microsoft-azure-batch"
    static public final String DEFAULT_OFFER = "centos-container"
    static public final String DEFAULT_VM_TYPE = "Standard_A3"
    static public final OSType DEFAULT_OS_TYPE = OSType.LINUX

    String publisher
    String offer
    OSType osType = DEFAULT_OS_TYPE
    VerificationType verification = VerificationType.VERIFIED

    String vmType
    Integer vmCount = 1
    boolean autoScale
    String scaleFormula

    String schedulePolicy // spread | pack

    AzPoolOpts() {
        this(Collections.emptyMap())
    }

    AzPoolOpts(Map opts) {
        this.publisher = opts.publisher ?: DEFAULT_PUBLISHER
        this.offer = opts.offer ?: DEFAULT_OFFER
        this.vmType = opts.vmType ?: DEFAULT_VM_TYPE
        this.vmCount = opts.vmCount as Integer ?: 1
        this.autoScale = opts.autoScale as boolean
        this.scaleFormula = opts.scaleFormula
        this.schedulePolicy = opts.schedulePolicy
    }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(publisher)
        hasher.putUnencodedChars(offer)
        hasher.putUnencodedChars(vmType)
        hasher.putInt(vmCount)
        hasher.putBoolean(autoScale)
        hasher.putUnencodedChars(scaleFormula ?: '')
        hasher.putUnencodedChars(schedulePolicy ?: '')
        return hasher
    }
}
