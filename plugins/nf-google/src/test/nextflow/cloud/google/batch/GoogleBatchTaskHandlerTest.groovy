/*
 * Copyright 2022, Google Inc.
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
 *
 */

package nextflow.cloud.google.batch

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.google.batch.model.ProvisioningModel
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.util.CpuUnit
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleBatchTaskHandlerTest extends Specification {

    def 'should create network policy' () {
        given:
        def handler = new GoogleBatchTaskHandler()

        when:
        def policy = handler.networkPolicy(Mock(BatchConfig))
        then:
        policy == null

        when:
        def config2 = Mock(BatchConfig) {getNetwork() >> 'net-abc'; getSubnetwork() >> 'sub-123' }
        def policy2 = handler.networkPolicy(config2)
        then:
        policy2.getNetworkInterfaces().size() == 1
        and:
        policy2.getNetworkInterfaces().get(0).getNetwork() == 'net-abc'
        policy2.getNetworkInterfaces().get(0).getSubnetwork() == 'sub-123'
        !policy2.getNetworkInterfaces().get(0).noExternalIpAddress

        when:
        def config3 = Mock(BatchConfig) {getUsePrivateAddress() >> true }
        def policy3 = handler.networkPolicy(config3)
        then:
        policy3.getNetworkInterfaces().size() == 1
        and:
        policy3.getNetworkInterfaces().get(0).noExternalIpAddress
    }

    def 'should create submit request' () {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        def CONTAINER_IMAGE = 'debian:latest'
        def CONTAINER_OPTS = '--foo'
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig)
        }
        and:
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        def task = Mock(TaskRun) {
            toTaskBean() >> bean
            getHashLog() >> 'abcd1234'
            getWorkDir() >> WORK_DIR
            getContainer() >> CONTAINER_IMAGE
            getConfig() >> Mock(TaskConfig) {
                getCpuUnits() >> CpuUnit.of(2)
                getContainerOptions() >> CONTAINER_OPTS
            }
        }

        and:
        def handler = new GoogleBatchTaskHandler(task, exec)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        req.getTaskGroups().size() == 1
        and:
        def group = req.getTaskGroups().get(0)
        and:
        group.taskSpec.computeResource.cpuMilli == 2_000
        group.taskSpec.computeResource.memoryMib == null
        group.taskSpec.maxRunDuration == null
        and:
        group.taskSpec.runnables.size() == 1
        def runnable = group.taskSpec.runnables.get(0)
        runnable.container.getCommands().join(' ') == '/bin/bash -o pipefail -c trap "{ cp .command.log /mnt/foo/scratch/.command.log; }" ERR; /bin/bash /mnt/foo/scratch/.command.run 2>&1 | tee .command.log'
        runnable.container.getImageUri() == CONTAINER_IMAGE
        runnable.container.getOptions() == CONTAINER_OPTS
        runnable.container.getVolumes() == ['/mnt/foo/scratch:/mnt/foo/scratch:rw']
        and:
        !req.allocationPolicy.network
        !req.allocationPolicy.instancePolicy().machineType
    }

    def 'should create submit request/2' () {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        and:
        def CONTAINER_IMAGE = 'ubuntu:22.1'
        def CONTAINER_OPTS = '--this --that'
        def CPUS = CpuUnit.of(4)
        def MEM = MemoryUnit.of('8 GB')
        def TIMEOUT = Duration.of('1 hour')
        def MACHINE_TYPE = 'vm-type-2'
        and:
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig) {
                getSpot() >> true
                getNetwork() >> 'net-1'
                getSubnetwork() >> 'subnet-1'
                getUsePrivateAddress() >> true
            }
        }
        and:
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        def task = Mock(TaskRun) {
            toTaskBean() >> bean
            getHashLog() >> 'abcd1234'
            getWorkDir() >> WORK_DIR
            getContainer() >> CONTAINER_IMAGE
            getConfig() >> Mock(TaskConfig) {
                getCpuUnits() >> CPUS
                getMemory() >> MEM
                getTime() >> TIMEOUT
                getMachineType() >> MACHINE_TYPE
                getContainerOptions() >> CONTAINER_OPTS
            }
        }

        and:
        def handler = new GoogleBatchTaskHandler(task, exec)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        req.getTaskGroups().size() == 1
        and:
        def group = req.getTaskGroups().get(0)
        and:
        group.taskSpec.computeResource.cpuMilli == CPUS.toMillis()
        group.taskSpec.computeResource.memoryMib == MEM.toMega()
        group.taskSpec.maxRunDuration == TIMEOUT.seconds + 's'
        and:
        group.taskSpec.runnables.size() == 1
        def runnable = group.taskSpec.runnables.get(0)
        runnable.container.getCommands().join(' ') == '/bin/bash -o pipefail -c trap "{ cp .command.log /mnt/foo/scratch/.command.log; }" ERR; /bin/bash /mnt/foo/scratch/.command.run 2>&1 | tee .command.log'
        runnable.container.getImageUri() == CONTAINER_IMAGE
        runnable.container.getOptions() == CONTAINER_OPTS
        runnable.container.getVolumes() == ['/mnt/foo/scratch:/mnt/foo/scratch:rw']
        and:
        req.allocationPolicy.instancePolicy().provisioningModel == ProvisioningModel.SPOT
        req.allocationPolicy.instancePolicy().machineType == MACHINE_TYPE
        and:
        def netInterface = req.allocationPolicy.network.getNetworkInterfaces().get(0)
        netInterface.network == 'net-1'
        netInterface.subnetwork == 'subnet-1'
        netInterface.noExternalIpAddress
    }
}
