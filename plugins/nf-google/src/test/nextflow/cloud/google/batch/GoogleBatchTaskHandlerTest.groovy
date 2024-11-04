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

import java.nio.file.Path

import com.google.cloud.batch.v1.GCS
import com.google.cloud.batch.v1.StatusEvent
import com.google.cloud.batch.v1.TaskStatus
import com.google.cloud.batch.v1.Volume
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.cloud.types.PriceModel
import nextflow.executor.Executor
import nextflow.executor.res.AcceleratorResource
import nextflow.executor.res.DiskResource
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
        def GCS_VOL = Volume.newBuilder().setGcs(GCS.newBuilder().setRemotePath('foo').build() ).build()
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
        def mounts = ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        def volumes = [GCS_VOL]
        def launcher = new GoogleBatchLauncherSpecMock('bash .command.run', mounts, volumes)
        
        and:
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))

        when:
        def req = handler.newSubmitRequest(task, launcher)
        then:
        handler.fusionEnabled() >> false
        handler.findBestMachineType(_, false) >> null

        and:
        def taskGroup = req.getTaskGroups(0)
        def runnable = taskGroup.getTaskSpec().getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicyOrTemplate = allocationPolicy.getInstances(0)
        def instancePolicy = allocationPolicy.getInstances(0).getPolicy()
        and:
        taskGroup.getTaskSpec().getComputeResource().getBootDiskMib() == 0
        taskGroup.getTaskSpec().getComputeResource().getCpuMilli() == 2_000
        taskGroup.getTaskSpec().getComputeResource().getMemoryMib() == 0
        taskGroup.getTaskSpec().getMaxRunDuration().getSeconds() == 0
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c bash .command.run'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        !runnable.getContainer().getOptions()
        runnable.getContainer().getVolumesList() == ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        and:
        !instancePolicyOrTemplate.getInstanceTemplate()
        and:
        instancePolicy.getAcceleratorsCount() == 0
        instancePolicy.getDisksCount() == 0
        !instancePolicy.getMachineType()
        !instancePolicy.getMinCpuPlatform()
        instancePolicy.getProvisioningModel().toString() == 'PROVISIONING_MODEL_UNSPECIFIED'
        !instancePolicy.getBootDisk().getImage()
        and:
        allocationPolicy.getLocation().getAllowedLocationsCount() == 0
        allocationPolicy.getNetwork().getNetworkInterfacesCount() == 0
        and:
        req.getLogsPolicy().getDestination().toString() == 'CLOUD_LOGGING'
        and:
        taskGroup.getTaskSpec().getVolumesList().size()==1
        taskGroup.getTaskSpec().getVolumes(0) == GCS_VOL
    }

    def 'should create submit request with maximal spec' () {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        and:
        def ACCELERATOR = new AcceleratorResource(request: 1, type: 'nvidia-tesla-v100')
        def BOOT_DISK = MemoryUnit.of('10 GB')
        def BOOT_IMAGE = 'batch-debian'
        def CONTAINER_IMAGE = 'ubuntu:22.1'
        def CONTAINER_OPTS = '--this --that'
        def CPU_PLATFORM = 'Intel Skylake'
        def CPUS = 4
        def DISK = new DiskResource(request: '100 GB', type: 'pd-standard')
        def MACHINE_TYPE = 'vm-type-2'
        def MEM = MemoryUnit.of('8 GB')
        def TIMEOUT = Duration.of('1 hour')
        and:
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig) {
                getAllowedLocations() >> ['zones/us-central1-a', 'zones/us-central1-c']
                getBootDiskSize() >> BOOT_DISK
                getBootDiskImage() >> BOOT_IMAGE
                getCpuPlatform() >> CPU_PLATFORM
                getMaxSpotAttempts() >> 5
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
                getDiskResource() >> DISK
                getMachineType() >> MACHINE_TYPE
                getMemory() >> MEM
                getTime() >> TIMEOUT
                getResourceLabels() >> [foo: 'bar']
            }
        }
        and:
        def mounts = ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        def launcher = new GoogleBatchLauncherSpecMock( 'bash .command.run', mounts )

        and:
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))

        when:
        def req = handler.newSubmitRequest(task, launcher)
        then:
        handler.fusionEnabled() >> false
        handler.findBestMachineType(_, false) >> new GoogleBatchMachineTypeSelector.MachineType(type: MACHINE_TYPE, location: "location", priceModel: PriceModel.spot)

        and:
        def taskGroup = req.getTaskGroups(0)
        def taskSpec = taskGroup.getTaskSpec()
        def runnable = taskSpec.getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicy = allocationPolicy.getInstances(0).getPolicy()
        def networkInterface = allocationPolicy.getNetwork().getNetworkInterfaces(0)
        and:
        taskSpec.getComputeResource().getBootDiskMib() == BOOT_DISK.toMega()
        taskSpec.getComputeResource().getCpuMilli() == CPUS * 1_000
        taskSpec.getComputeResource().getMemoryMib() == MEM.toMega()
        taskSpec.getMaxRunDuration().getSeconds() == TIMEOUT.seconds
        taskSpec.getVolumes(0).getMountPath() == '/tmp'
        taskSpec.getMaxRetryCount() == 5
        taskSpec.getLifecyclePolicies(0).getActionCondition().getExitCodes(0) == 50001
        taskSpec.getLifecyclePolicies(0).getAction().toString() == 'RETRY_TASK'
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c bash .command.run'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        runnable.getContainer().getOptions() == '--this --that --privileged'
        runnable.getContainer().getVolumesList() == [
            '/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw',
            '/var/lib/nvidia/lib64:/usr/local/nvidia/lib64',
            '/var/lib/nvidia/bin:/usr/local/nvidia/bin'
        ]
        and:
        runnable.getEnvironment().getVariablesMap() == [:]
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
        instancePolicy.getBootDisk().getImage() == BOOT_IMAGE
        instancePolicy.getDisks(0).getNewDisk().getSizeGb() == DISK.request.toGiga()
        instancePolicy.getDisks(0).getNewDisk().getType() == DISK.type
        instancePolicy.getMachineType() == MACHINE_TYPE
        instancePolicy.getMinCpuPlatform() == CPU_PLATFORM
        instancePolicy.getProvisioningModel().toString() == 'SPOT'
        and:
        networkInterface.getNetwork() == 'net-1'
        networkInterface.getSubnetwork() == 'subnet-1'
        networkInterface.getNoExternalIpAddress() == true
        and:
        req.getLogsPolicy().getDestination().toString() == 'CLOUD_LOGGING'
        and:
        req.getLabelsMap() == [foo: 'bar']


        // with custom disk type
        when:
        req = handler.newSubmitRequest(task, launcher)
        then:
        task.getConfig().getDiskResource() >> new DiskResource(request: '100 GB')
        handler.fusionEnabled() >> false
        handler.findBestMachineType(_, false) >> null
        and:
        req.getTaskGroups(0).getTaskSpec().getComputeResource().getBootDiskMib() == 100 * 1024
    }

    def 'should use custom job name'() {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        def CONTAINER_IMAGE = 'debian:latest'
        and:
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        def sess = Mock(Session)
        def exec = new GoogleBatchExecutor(session: sess)

        and:
        def task = Mock(TaskRun) {
            toTaskBean() >> bean
            getHashLog() >> 'abcd1234'
            getWorkDir() >> WORK_DIR
            getContainer() >> CONTAINER_IMAGE
            getConfig() >> Mock(TaskConfig)
        }
        and:
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))

        when:
        def result = handler.customJobName(task)
        then:
        1 * sess.getExecConfigProp(_,'jobName', null) >> null
        and:
        result == null

        when:
        result = handler.customJobName(task)
        then:
        1 * sess.getExecConfigProp(_,'jobName', null) >> { return { "foo-${task.hashLog}" }  }
        and:
        result == 'foo-abcd1234'

    }

    def 'should use instance template' () {
        given:
        def GCS_VOL = Volume.newBuilder().setGcs(GCS.newBuilder().setRemotePath('foo').build() ).build()
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        def CONTAINER_IMAGE = 'debian:latest'
        def INSTANCE_TEMPLATE = 'instance-template'
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig) {
                getInstallGpuDrivers() >> true
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
                getCpus() >> 2
                getMachineType() >> "template://${INSTANCE_TEMPLATE}"
                getResourceLabels() >> [:]
            }
        }
        and:
        def mounts = ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        def volumes = [GCS_VOL]
        def launcher = new GoogleBatchLauncherSpecMock('bash .command.run', mounts, volumes)

        and:
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))

        when:
        def req = handler.newSubmitRequest(task, launcher)
        then:
        handler.fusionEnabled() >> false
        handler.findBestMachineType(_, false) >> null

        and:
        def taskGroup = req.getTaskGroups(0)
        def runnable = taskGroup.getTaskSpec().getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicyOrTemplate = allocationPolicy.getInstances(0)
        and:
        taskGroup.getTaskSpec().getComputeResource().getBootDiskMib() == 0
        taskGroup.getTaskSpec().getComputeResource().getCpuMilli() == 2_000
        taskGroup.getTaskSpec().getComputeResource().getMemoryMib() == 0
        taskGroup.getTaskSpec().getMaxRunDuration().getSeconds() == 0
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c bash .command.run'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        !runnable.getContainer().getOptions()
        runnable.getContainer().getVolumesList() == ['/mnt/disks/foo/scratch:/mnt/disks/foo/scratch:rw']
        and:
        instancePolicyOrTemplate.getInstallGpuDrivers() == true
        instancePolicyOrTemplate.getInstanceTemplate() == INSTANCE_TEMPLATE
        and:
        taskGroup.getTaskSpec().getVolumesCount() == 1
        taskGroup.getTaskSpec().getVolumes(0) == GCS_VOL
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
        handler.@taskId = '0'
        handler.@uid = '789'

        when:
        def trace = handler.getTraceRecord()
        then:
        handler.isCompleted() >> false
        and:
        trace.native_id == 'xyz-123/0/789'
        trace.executorName == 'google-batch'
    }

    def 'should create submit request with fusion enabled' () {
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
        def env = [FUSION_WORK: '/xyz']
        def launcher = new GoogleBatchLauncherSpecMock('bash .command.run', [], [], env)

        and:
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))

        when:
        def req = handler.newSubmitRequest(task, launcher)
        then:
        handler.fusionEnabled() >> true
        handler.findBestMachineType(_, true) >> null
        and:
        def taskGroup = req.getTaskGroups(0)
        def runnable = taskGroup.getTaskSpec().getRunnables(0)
        def allocationPolicy = req.getAllocationPolicy()
        def instancePolicy = allocationPolicy.getInstances(0).getPolicy()
        and:
        taskGroup.getTaskSpec().getComputeResource().getCpuMilli() == 2_000
        taskGroup.getTaskSpec().getVolumes(0).getMountPath() == '/tmp'
        and:
        runnable.getContainer().getCommandsList().join(' ') == '/bin/bash -o pipefail -c bash .command.run'
        runnable.getContainer().getImageUri() == CONTAINER_IMAGE
        runnable.getContainer().getOptions() == '--privileged'
        runnable.getContainer().getVolumesCount() == 0
        and:
        runnable.getEnvironment().getVariablesMap() == env
        and:
        instancePolicy.getAcceleratorsCount() == 0
        instancePolicy.getDisks(0).getNewDisk().getSizeGb() == 375
        instancePolicy.getDisks(0).getNewDisk().getType() == 'local-ssd'
        !instancePolicy.getMachineType()
        !instancePolicy.getMinCpuPlatform()
        instancePolicy.getProvisioningModel().toString() == 'PROVISIONING_MODEL_UNSPECIFIED'
        and:
        allocationPolicy.getLocation().getAllowedLocationsCount() == 0
        allocationPolicy.getNetwork().getNetworkInterfacesCount() == 0
        and:
        req.getLogsPolicy().getDestination().toString() == 'CLOUD_LOGGING'
    }

    def 'should not set wildcard expressions as machine type'() {
        given:
        def WORK_DIR = CloudStorageFileSystem.forBucket('foo').getPath('/scratch')
        def CONTAINER_IMAGE = 'debian:latest'
        def exec = Mock(GoogleBatchExecutor) {
            getConfig() >> Mock(BatchConfig)
        }
        def bean = new TaskBean(workDir: WORK_DIR, inputFiles: [:])
        def task = Mock(TaskRun) {
            toTaskBean() >> bean
            getHashLog() >> 'abcd1234'
            getWorkDir() >> WORK_DIR
            getContainer() >> CONTAINER_IMAGE
            getConfig() >> Mock(TaskConfig) {
                getCpus() >> 2
                getResourceLabels() >> [:]
                getMachineType() >> "n1-*,n2-*"
            }
        }
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))
        def env = [FUSION_WORK: '/xyz']
        def launcher = new GoogleBatchLauncherSpecMock('bash .command.run', [], [], env)

        when:
        def req = handler.newSubmitRequest(task, launcher)
        then:
        handler.fusionEnabled() >> false
        handler.findBestMachineType(_, _) >> null
        and:
        req.getAllocationPolicy().getInstances(0).policy.getMachineType() == ""

    }

    TaskStatus makeTaskStatus(String desc) {
        TaskStatus.newBuilder()
            .addStatusEvents(
                StatusEvent.newBuilder()
                    .setDescription(desc)
            )
            .build()
    }

    def 'should detect spot failures from status event'() {
        given:
        def jobId = 'job-id'
        def taskId = 'task-id'
        def client = Mock(BatchClient)
        def task = Mock(TaskRun) {
            lazyName() >> 'foo (1)'
        }
        def handler = Spy(new GoogleBatchTaskHandler(jobId: jobId, taskId: taskId, client: client, task: task))

        when:
        client.getTaskStatus(jobId, taskId) >>> [
            makeTaskStatus('Task failed due to Spot VM preemption with exit code 50001.'),
            makeTaskStatus('Task succeeded')
        ]
        then:
        handler.getJobError().message == "Task failed due to Spot VM preemption with exit code 50001."
        handler.getJobError() == null
    }

    def 'should find best instance type' () {
        given:
        def workDir = Path.of('/work/dir')
        def client = Mock(BatchClient)
        def task = Mock(TaskRun) {
             hashLog >> '1234567890'
             getWorkDir() >> workDir
        }
        def exec = Mock(GoogleBatchExecutor) {
            getClient() >> client
            getConfig() >> Mock(BatchConfig) { getSpot()>>false }
            isCloudinfoEnabled() >> true
        }
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))
        and:
        def config = Mock(TaskConfig)
        def machineType = GroovyMock(GoogleBatchMachineTypeSelector.MachineType)

        when:
        def result = handler.findBestMachineType(config, false)
        then:
        1 * handler.bestMachineType0(_,_,_,_,_,_) >> machineType
        and:
        result == machineType
    }

    def 'should disable cloudinfo' () {
        given:
        SysEnv.push(NXF_CLOUDINFO_ENABLED: 'false')

        def workDir = Path.of('/work/dir')
        def client = Mock(BatchClient)
        def task = Mock(TaskRun) {
            hashLog >> '1234567890'
            getWorkDir() >> workDir
        }
        def exec = Mock(GoogleBatchExecutor) {
            getClient() >> client
            getConfig() >> Mock(BatchConfig) { getSpot()>>false }
            isCloudinfoEnabled() >> false
        }
        def handler = Spy(new GoogleBatchTaskHandler(task, exec))
        and:
        def config = Mock(TaskConfig)

        when:
        def result = handler.findBestMachineType(config, false)
        then:
        0 * handler.bestMachineType0(_,_,_,_,_,_) >> null
        and:
        result == null

        cleanup:
        SysEnv.pop()
    }

    def 'should kill a job' () {
        given:
        def client = Mock(BatchClient)
        def executor = Mock(GoogleBatchExecutor)
        def task = Mock(TaskRun)
        def handler = Spy(GoogleBatchTaskHandler)
        handler.@executor = executor
        handler.@client = client
        handler.task = task

        when:
        handler.@jobId = 'job1'
        handler.kill()
        then:
        handler.isActive() >> false
        0 * executor.shouldDeleteJob('job1') >> true
        and:
        0 * client.deleteJob('job1') >> null

        when:
        handler.@jobId = 'job1'
        handler.kill()
        then:
        handler.isActive() >> true
        1 * executor.shouldDeleteJob('job1') >> true
        and:
        1 * client.deleteJob('job1') >> null

        when:
        handler.@jobId = 'job1'
        handler.kill()
        then:
        handler.isActive() >> true
        1 * executor.shouldDeleteJob('job1') >> false
        and:
        0 * client.deleteJob('job1') >> null
    }
}
