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

package nextflow.cloud.azure.batch

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.transform.builder.Builder
import nextflow.cloud.azure.config.AzPoolOpts

/**
 * Model the spec of Azure VM pool
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Builder
@ToString(includeNames = true, includePackage = false)
@EqualsAndHashCode
@CompileStatic
class AzVmPoolSpec {
    String poolId
    AzVmType vmType
    AzPoolOpts opts
}
