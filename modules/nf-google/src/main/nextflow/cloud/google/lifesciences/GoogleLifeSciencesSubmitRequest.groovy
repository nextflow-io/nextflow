/*
 * Copyright 2019, Google Inc
 * Copyright 2018, WuxiNextcode
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

package nextflow.cloud.google.lifesciences

import java.nio.file.Path

import com.google.api.services.lifesciences.v2beta.model.Mount
import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.executor.res.AcceleratorResource
/**
 * Models Google pipeline request for a Nextflow task executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
class GoogleLifeSciencesSubmitRequest {

    String machineType

    String project

    List<String> zone

    List<String> region

    String diskName

    Integer diskSizeGb

    boolean preemptible

    String taskName

    String containerImage

    Mount sharedMount

    AcceleratorResource accelerator

    String location

    Path workDir

    Integer bootDiskSizeGb
    
    String cpuPlatform

    String entryPoint

    boolean usePrivateAddress
}
