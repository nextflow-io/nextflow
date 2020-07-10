/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.nio.file.Paths

import com.google.api.services.lifesciences.v2beta.model.Event
import com.google.api.services.lifesciences.v2beta.model.Metadata
import com.google.api.services.lifesciences.v2beta.model.Operation
import nextflow.Session
import nextflow.cloud.google.GoogleSpecification
import nextflow.cloud.types.PriceModel
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.script.BaseScript
import nextflow.script.ProcessConfig
import nextflow.util.CacheHelper
import nextflow.util.MemoryUnit
import spock.lang.Shared

class GoogleLifeSciencesTaskHandlerTest extends GoogleSpecification {

    @Shared
    GoogleLifeSciencesConfig pipeConfig = new GoogleLifeSciencesConfig(
            project: "testProject",
            zones: ["testZone"],
            regions: ["testRegion"])

    @Shared
    UUID uuid = new UUID(4,4)

    @Shared
    Session stubSession = Stub {
        getUniqueId() >> uuid
        getConfig() >> ["cloud" : ["instanceType": "instanceType"]]
    }

    GoogleLifeSciencesExecutor stubExecutor = GroovyStub() {
        getSession() >> stubSession
        getHelper() >> GroovyMock(GoogleLifeSciencesHelper)
    }

    TaskRun stubTaskRunner = GroovyStub {
        getName() >> "testName"
        getId() >> new TaskId(12345)
        getContainer() >> "testContainer"
        getScript() >> "echo testScript"
        getWorkDir() >> File.createTempDir().toPath()
        getHash() >> { CacheHelper.hasher('dummy').hash() }
    }

    def 'should throw an error if container is not specified'() {
        given:
        def noContainerTaskRunner = Stub(TaskRun) {
            getName() >> "noContainer"
        }

        when: 'handler is constructed'
        def handler = new GoogleLifeSciencesTaskHandler(noContainerTaskRunner,stubExecutor)

        then: 'we should get an error stating that container definition is missing'
        def error = thrown(ProcessUnrecoverableException)
        error.getMessage() == "No container image specified for process $noContainerTaskRunner.name -- Either specify the container to use in the process definition or with 'process.container' value in your config"
        !handler

    }

    def 'should construct correctly'() {
        when:
        def handler = new GoogleLifeSciencesTaskHandler(stubTaskRunner, stubExecutor)

        then:
        handler.task.container == "testContainer"
    }



    def 'should submit a task'() {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'foo'

        def handler = Spy(GoogleLifeSciencesTaskHandler)
        handler.task = task

        def req = Mock(GoogleLifeSciencesSubmitRequest)
        def operation = new Operation()

        when:
        handler.submit()

        then:
        1 * handler.createTaskWrapper() >> null
        1 * handler.createPipelineRequest() >> req
        1 * handler.submitPipeline(req) >> operation
        1 * handler.getPipelineIdFromOp(operation) >> '123'

        handler.operation == operation
        handler.pipelineId == '123'
        handler.status == TaskStatus.SUBMITTED

    }

    def 'should check submitPipeline' () {
        given:
        def helper = Mock(GoogleLifeSciencesHelper)
        def executor = new GoogleLifeSciencesExecutor(helper: helper)
        def handler = new GoogleLifeSciencesTaskHandler(executor: executor, helper: helper)

        def operation = new Operation()
        def request = Mock(GoogleLifeSciencesSubmitRequest)

        when:
        def op = handler.submitPipeline(request)
        then:
        1 * helper.submitPipeline(request) >> operation
        op == operation
    }

    def 'should create pipeline request' () {
        given:
        def workDir = mockGsPath('gs://my-bucket/work/dir')
        and:
        def executor = Mock(GoogleLifeSciencesExecutor) {
            getHelper() >> Mock(GoogleLifeSciencesHelper)
            getConfig() >> {
                Mock(GoogleLifeSciencesConfig) {
                    getProject() >> 'my-project'
                    getZones() >> ['my-zone']
                    getRegions() >> ['my-region']
                }
            }
        }
        and:
        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'

        def handler = new GoogleLifeSciencesTaskHandler(
                executor: executor,
                task: task )

        when:
        def req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(machineType: 'n1-1234')
        and:
        req.machineType == 'n1-1234'
        req.project == 'my-project'
        req.zone == ['my-zone']
        req.region == ['my-region']
        req.diskName == GoogleLifeSciencesTaskHandler.DEFAULT_DISK_NAME
        req.diskSizeGb == null
        !req.preemptible
        req.taskName == "nf-bad893071e9130b866d43a4fcabb95b6"
        req.containerImage == 'my/image'
        req.workDir.toUriString() == 'gs://my-bucket/work/dir'
        req.sharedMount.getPath() == '/work/dir'
        req.sharedMount.getDisk() == GoogleLifeSciencesTaskHandler.DEFAULT_DISK_NAME
        !req.sharedMount.getReadOnly()
        req.bootDiskSizeGb == null
        req.cpuPlatform == null
        req.entryPoint == GoogleLifeSciencesConfig.DEFAULT_ENTRY_POINT
        !req.usePrivateAddress

        when:
        req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(containerOptions: [entrypoint:'/bin/foo'])
        and:
        req.entryPoint == '/bin/foo'

        when:
        req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(containerOptions: [entrypoint:null])
        and:
        req.entryPoint == null

    }

    def 'should create pipeline request/2' () {
        given:
        def workDir = mockGsPath('gs://my-bucket/work/dir')
        and:
        def executor = Mock(GoogleLifeSciencesExecutor) {
            getHelper() >> Mock(GoogleLifeSciencesHelper)
            getConfig() >> {
                Mock(GoogleLifeSciencesConfig) {
                    getProject() >> 'my-project'
                    getZones() >> ['my-zone']
                    getRegions() >> ['my-region']
                    getPreemptible() >> true
                    getBootDiskSize() >> MemoryUnit.of('20 GB')
                    getUsePrivateAddress() >> true
                    getCpuPlatform() >> 'Intel Skylake'
                }
            }
        }
        and:
        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'

        def handler = new GoogleLifeSciencesTaskHandler(
                executor: executor,
                task: task )

        when:
        def req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(disk: '250 GB', machineType: 'n1-1234')
        and:
        req.machineType == 'n1-1234'
        req.project == 'my-project'
        req.zone == ['my-zone']
        req.region == ['my-region']
        req.diskName == GoogleLifeSciencesTaskHandler.DEFAULT_DISK_NAME
        req.diskSizeGb == 250
        req.preemptible
        req.taskName == "nf-bad893071e9130b866d43a4fcabb95b6"
        req.containerImage == 'my/image'
        req.workDir.toUriString() == 'gs://my-bucket/work/dir'
        req.sharedMount.getPath() == '/work/dir'
        req.sharedMount.getDisk() == GoogleLifeSciencesTaskHandler.DEFAULT_DISK_NAME
        !req.sharedMount.getReadOnly()
        req.bootDiskSizeGb == 20
        req.cpuPlatform =='Intel Skylake'
        req.entryPoint == GoogleLifeSciencesConfig.DEFAULT_ENTRY_POINT
        req.usePrivateAddress

        when:
        req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(containerOptions: [entrypoint:'/bin/foo'])
        and:
        req.entryPoint == '/bin/foo'

        when:
        req = handler.createPipelineRequest()
        then:
        task.getConfig() >> new TaskConfig(containerOptions: [entrypoint:null])
        and:
        req.entryPoint == null

    }

    
    def 'should check if it is running'(){
        given:
        // -- task
        def task = Mock(TaskRun)
        task.name >> 'google-task'
        task.getWorkDir() >> Paths.get('/work/dir')

        // -- executor
        def helper = Mock(GoogleLifeSciencesHelper)
        def executor = new GoogleLifeSciencesExecutor(helper: helper)

        def operation = new Operation()

        // -- handler
        def handler = Spy(GoogleLifeSciencesTaskHandler)
        handler.executor = executor
        handler.task = task
        handler.operation = operation
        handler.helper = helper

        when:
        def result = handler.checkIfRunning()
        then:
        1 * handler.isSubmitted() >> false
        0 * handler.executor.helper.checkOperationStatus(_)
        handler.status == TaskStatus.NEW
        result == false

        when:
        result = handler.checkIfRunning()
        then:
        1 * handler.isSubmitted() >> true
        1 * handler.executor.helper.checkOperationStatus(_) >> { new Operation() }
        handler.status == TaskStatus.RUNNING
        result == true

        when:
        result = handler.checkIfRunning()
        then:
        1 * handler.isSubmitted() >> false
        0 * handler.executor.helper.checkOperationStatus(_)
        handler.status == TaskStatus.RUNNING
        result == false

    }


    def 'should check if it is complete'() {
        given:
        // -- task
        def task = Mock(TaskRun)
        task.name >> 'google-task'
        task.getWorkDir() >> Paths.get('/work/dir')

        // -- executor
        def helper = Mock(GoogleLifeSciencesHelper)
        def executor = new GoogleLifeSciencesExecutor(helper:helper)

        // -- handler
        def handler = Spy(GoogleLifeSciencesTaskHandler)
        handler.executor = executor
        handler.task = task
        handler.helper = helper

        when:
        def isComplete = handler.checkIfCompleted()
        then:
        1 * handler.isRunning() >> false
        0 * helper.checkOperationStatus(_)
        !isComplete

        when:
        isComplete = handler.checkIfCompleted()
        then:
        1 * handler.isRunning() >> true
        1 * helper.checkOperationStatus(_)  >> { new Operation().setName("incomplete").setDone(false) }
        !isComplete

        when:
        isComplete = handler.checkIfCompleted()
        then:
        1 * handler.readExitFile() >> 0
        1 * handler.isRunning() >> true
        1 * helper.checkOperationStatus(_) >> { new Operation().setName("complete").setDone(true) }
        isComplete
    }

    def 'should get jobId from operation' () {
        given:
        def operation = new Operation().setName('projects/rare-lattice-222412/operations/16737869387120678662')
        def handler = [:] as GoogleLifeSciencesTaskHandler
        expect:
        handler.getPipelineIdFromOp(operation) == '16737869387120678662'
    }

    def 'should get events from operation' () {
        given:
        def handler = [:] as GoogleLifeSciencesTaskHandler

        when:
        def op = new Operation()
        def events = handler.getEventsFromOp(op)
        then:
        events == []

        when:
        def e1 = new Event().setDescription('foo').setTimestamp('2018-12-15T12:50:30.743109Z')
        op.setMetadata(new Metadata().setEvents([e1]))
        events = handler.getEventsFromOp(op)
        then:
        events == [e1]

        when:
        def e2 = new Event().setDescription('bar').setTimestamp('2018-12-15T12:52:30.743109Z')
        def e3 = new Event().setDescription('baz').setTimestamp('2018-12-15T12:55:30.743109Z')
        op.setMetadata(new Metadata().setEvents([e3, e2, e1]))
        events = handler.getEventsFromOp(op)
        then:
        events == [e2, e3]

    }

    def 'should create disk mount'() {
        given:
        def diskName = "testDisk"
        def mountPath = "testPath"
        def readOnly = true
        and:
        def handler = Spy(GoogleLifeSciencesTaskHandler)

        when:
        def mount1 = handler.configureMount(diskName,mountPath,readOnly)
        then:
        with(mount1) {
            getDisk() == diskName
            getPath() == mountPath
            getReadOnly() == readOnly
        }

        when:
        def mount2 = handler.configureMount(diskName,mountPath)
        then:
        with(mount2) {
            getDisk() == diskName
            getPath() == mountPath
            !getReadOnly()
        }

    }

    def 'should create trace record'() {
        given:
        def executor = Mock(GoogleLifeSciencesExecutor)  {
            getConfig() >> new GoogleLifeSciencesConfig(location: 'eu-east-1', preemptible: true)
            getName() >> 'google-lifesciences'
        }
        and:
        def processor = Mock(TaskProcessor)
        processor.getExecutor() >> executor
        processor.getName() >> 'foo'
        processor.getConfig() >> new ProcessConfig(Mock(BaseScript))
        and:
        def task = Mock(TaskRun)
        task.getProcessor() >> processor
        task.getConfig() >> Mock(TaskConfig) { getMachineType() >> 'm1.large' }
        and:
        def handler = Spy(GoogleLifeSciencesTaskHandler)
        handler.task = task
        handler.pipelineId = 'xyz-123'
        handler.executor = executor
        handler.assignedZone = 'eu-east-1'

        when:
        def record = handler.getTraceRecord()
        then:
        record.get('native_id') == 'xyz-123'
        record.getExecutorName() == 'google-lifesciences'
        record.machineInfo.type == 'm1.large'
        record.machineInfo.zone == 'eu-east-1'
        record.machineInfo.priceModel == PriceModel.spot
    }

}
