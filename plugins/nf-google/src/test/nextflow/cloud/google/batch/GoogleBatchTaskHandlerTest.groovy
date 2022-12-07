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
import nextflow.executor.Executor
import nextflow.executor.res.AcceleratorResource
import nextflow.processor.TaskBean
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleBatchTaskHandlerTest extends Specification {

    def 'should create submit request with minimal spec' () {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        def CONTAINER_IMAGE = 'debian:latest'
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
                getCpus() >> 2
                getResourceLabels() >> [:]
            }
        }

        and:
        def handler = new GoogleBatchTaskHandler(task, exec)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        def taskGroup = req.getTaskGroups(0)
        def runnable = taskGroup.getTaskSpec().getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicy = allocationPolicy.getInstances(0).getPolicy()
        and:
        taskGroup.getTaskSpec().getComputeResource().getBootDiskMib() == 0
        taskGroup.getTaskSpec().getComputeResource().getCpuMilli() == 2_000
        taskGroup.getTaskSpec().getComputeResource().getMemoryMib() == 0
        taskGroup.getTaskSpec().getMaxRunDuration().getSeconds() == 0
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c trap "{ cp .command.log /mnt/disks/foo/scratch/.command.log; }" ERR; /bin/bash /mnt/disks/foo/scratch/.command.run 2>&1 | tee .command.log'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        runnable.getContainer().getOptions() == ''
        runnable.getContainer().getVolumesList() == ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        and:
        instancePolicy.getAcceleratorsCount() == 0
        instancePolicy.getMachineType() == ''
        instancePolicy.getMinCpuPlatform() == ''
        instancePolicy.getProvisioningModel().toString() == 'PROVISIONING_MODEL_UNSPECIFIED'
        and:
        allocationPolicy.getLocation().getAllowedLocationsCount() == 0
        allocationPolicy.getNetwork().getNetworkInterfacesCount() == 0
        and:
        req.getLogsPolicy().getDestination().toString() == 'CLOUD_LOGGING'
    }

    def 'should create submit request with maximal spec' () {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        and:
        def ACCELERATOR = new AcceleratorResource(request: 1, type: 'nvidia-tesla-v100')
        def BOOT_DISK = MemoryUnit.of('10 GB')
        def CONTAINER_IMAGE = 'ubuntu:22.1'
        def CONTAINER_OPTS = '--this --that'
        def CPU_PLATFORM = 'Intel Skylake'
        def CPUS = 4
        def DISK = MemoryUnit.of('50 GB')
        def MACHINE_TYPE = 'vm-type-2'
        def MEM = MemoryUnit.of('8 GB')
        def TIMEOUT = Duration.of('1 hour')
        and:
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig) {
                getAllowedLocations() >> ['zones/us-central1-a', 'zones/us-central1-c']
                getBootDiskSize() >> BOOT_DISK
                getCpuPlatform() >> CPU_PLATFORM
                getSpot() >> true
                getNetwork() >> 'net-1'
                getServiceAccountEmail() >> 'foo@bar.baz'
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
                getAccelerator() >> ACCELERATOR
                getContainerOptions() >> CONTAINER_OPTS
                getCpus() >> CPUS
                getDisk() >> DISK
                getMachineType() >> MACHINE_TYPE
                getMemory() >> MEM
                getTime() >> TIMEOUT
                getResourceLabels() >> [foo: 'bar']
            }
        }

        and:
        def handler = new GoogleBatchTaskHandler(task, exec)

        when:
        def req = handler.newSubmitRequest(task)
        then:
        def taskGroup = req.getTaskGroups(0)
        def runnable = taskGroup.getTaskSpec().getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicy = allocationPolicy.getInstances(0).getPolicy()
        def networkInterface = allocationPolicy.getNetwork().getNetworkInterfaces(0)
        and:
        taskGroup.getTaskSpec().getComputeResource().getBootDiskMib() == DISK.toMega()
        taskGroup.getTaskSpec().getComputeResource().getCpuMilli() == CPUS * 1_000
        taskGroup.getTaskSpec().getComputeResource().getMemoryMib() == MEM.toMega()
        taskGroup.getTaskSpec().getMaxRunDuration().getSeconds() == TIMEOUT.seconds
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c trap "{ cp .command.log /mnt/disks/foo/scratch/.command.log; }" ERR; /bin/bash /mnt/disks/foo/scratch/.command.run 2>&1 | tee .command.log'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        runnable.getContainer().getOptions() == '--this --that --privileged'
        runnable.getContainer().getVolumesList() == [
            '/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw',
            '/var/lib/nvidia/lib64:/usr/local/nvidia/lib64',
            '/var/lib/nvidia/bin:/usr/local/nvidia/bin'
        ]
        and:
        allocationPolicy.getLocation().getAllowedLocationsCount() == 2
        allocationPolicy.getLocation().getAllowedLocations(0) == 'zones/us-central1-a'
        allocationPolicy.getLocation().getAllowedLocations(1) == 'zones/us-central1-c'
        allocationPolicy.getInstances(0).getInstallGpuDrivers() == true
        allocationPolicy.getLabelsMap() == [foo: 'bar']
        allocationPolicy.getServiceAccount().getEmail() == 'foo@bar.baz'
        and:
        instancePolicy.getAccelerators(0).getCount() == 1
        instancePolicy.getAccelerators(0).getType() == ACCELERATOR.type
        instancePolicy.getMachineType() == MACHINE_TYPE
        instancePolicy.getMinCpuPlatform() == CPU_PLATFORM
        instancePolicy.getProvisioningModel().toString() == 'SPOT'
        and:
        networkInterface.getNetwork() == 'net-1'
        networkInterface.getSubnetwork() == 'subnet-1'
        networkInterface.getNoExternalIpAddress() == true
        and:
        req.getLogsPolicy().getDestination().toString() == 'CLOUD_LOGGING'
    }

    def 'should create the trace record' () {
        given:
        def exec = Mock(Executor) { getName() >> 'google-batch' }
        def processor = Mock(TaskProcessor) {
            getExecutor() >> exec
            getName() >> 'foo'
            getConfig() >> new ProcessConfig(Mock(BaseScript))
        }
        and:
        def task = Mock(TaskRun)
        task.getProcessor() >> processor
        task.getConfig() >> GroovyMock(TaskConfig)
        and:
        def handler = Spy(GoogleBatchTaskHandler)
        handler.task = task
        handler.@jobId = 'xyz-123'
        handler.@uid = '789'

        when:
        def trace = handler.getTraceRecord()
        then:
        handler.isCompleted() >> false
        and:
        trace.native_id == 'xyz-123/789'
        trace.executorName == 'google-batch'
    }
}
