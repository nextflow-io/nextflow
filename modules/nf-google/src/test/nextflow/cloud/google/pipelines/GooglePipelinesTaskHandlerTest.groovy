/*
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
package nextflow.cloud.google.pipelines

import spock.lang.Shared

import java.nio.file.Paths

import com.google.api.services.genomics.v2alpha1.model.Event
import com.google.api.services.genomics.v2alpha1.model.Metadata
import com.google.api.services.genomics.v2alpha1.model.Operation
import nextflow.Session
import nextflow.cloud.google.GoogleSpecification
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CacheHelper

class GooglePipelinesTaskHandlerTest extends GoogleSpecification {

    @Shared
    GooglePipelinesConfiguration pipeConfig = new GooglePipelinesConfiguration("testProject",["testZone"],["testRegion"])

    @Shared
    UUID uuid = new UUID(4,4)

    @Shared
    Session stubSession = Stub {
        getUniqueId() >> uuid
        getConfig() >> ["cloud" : ["instanceType": "instanceType"]]
    }

    GooglePipelinesExecutor stubExecutor = GroovyStub() {
        getSession() >> stubSession
        getHelper() >> GroovyMock(GooglePipelinesHelper)
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
        def handler = new GooglePipelinesTaskHandler(noContainerTaskRunner,stubExecutor,pipeConfig)

        then: 'we should get an error stating that container definition is missing'
        def error = thrown(ProcessUnrecoverableException)
        error.getMessage() == "No container image specified for process $noContainerTaskRunner.name -- Either specify the container to use in the process definition or with 'process.container' value in your config"
        !handler

    }

    def 'should construct correctly'() {
        when:
        def handler = new GooglePipelinesTaskHandler(stubTaskRunner, stubExecutor,pipeConfig)

        then:
        handler.task.container == "testContainer"
    }



    def 'should submit a task'() {
        given:
        def task = Mock(TaskRun)
        task.getName() >> 'foo'

        def handler = Spy(GooglePipelinesTaskHandler)
        handler.task = task

        def req = Mock(GooglePipelinesSubmitRequest)
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
        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)
        def handler = new GooglePipelinesTaskHandler(executor: executor)

        def operation = new Operation()
        def request = Mock(GooglePipelinesSubmitRequest)

        when:
        def op = handler.submitPipeline(request)
        then:
        1 * executor.helper.submitPipeline(request) >> operation
        op == operation
    }

    def 'should create pipeline request' () {
        given:

        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)
        def config = Mock(GooglePipelinesConfiguration)
        def workDir = mockGsPath('gs://my-bucket/work/dir')

        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'
        task.getConfig() >> new TaskConfig(disk: '250 GB')

        def handler = new GooglePipelinesTaskHandler(
                pipelineConfiguration: config,
                executor: executor,
                task: task,
                machineType: 'n1-1234'
        )

        when:
        def req = handler.createPipelineRequest()

        then:
        config.getProject() >> 'my-project'
        config.getZone() >> ['my-zone']
        config.getRegion() >> ['my-region']
        config.getPreemptible() >> true
        // chek request object
        req.machineType == 'n1-1234'
        req.project == 'my-project'
        req.zone == ['my-zone']
        req.region == ['my-region']
        req.diskName == GooglePipelinesTaskHandler.diskName
        req.diskSizeGb == 250
        req.preemptible == true
        req.taskName == "nf-bad893071e9130b866d43a4fcabb95b6"
        req.containerImage == 'my/image'
        req.fileCopyImage == GooglePipelinesTaskHandler.fileCopyImage
        // check staging script
        req.stagingScript == 'mkdir -p /work/dir'
        // check main script
        req.mainScript == 'cd /work/dir; bash .command.run 2>&1 | tee -a .command.log'
        // check unstaging script
        req.unstagingScript.tokenize(';')[0] == 'cd /work/dir'
        req.unstagingScript.tokenize(';')[1] == ' [[ $GOOGLE_PIPELINE_FAILED == 1 || $NXF_DEBUG ]] && gsutil -m -q cp -R /google/ gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[2] == ' [[ -f .command.trace ]] && gsutil -m -q cp -R .command.trace gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[3] == ' gsutil -m -q cp -R .command.err gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[4] == ' gsutil -m -q cp -R .command.out gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[5] == ' gsutil -m -q cp -R .command.log gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[6] == ' gsutil -m -q cp -R .exitcode gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';').size() == 7
    }

    def 'should create pipeline request with stage and unstage commands' () {
        given:

        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)
        def config = Mock(GooglePipelinesConfiguration)
        def workDir = mockGsPath('gs://my-bucket/work/dir')

        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'
        task.getConfig() >> new TaskConfig()

        def handler = new GooglePipelinesTaskHandler(
                pipelineConfiguration: config,
                executor: executor,
                machineType: 'n1-1234',
                task: task,
                stagingCommands: ['alpha', 'beta', 'delta'],
                unstagingCommands: ['foo', 'bar']
        )

        when:
        def req = handler.createPipelineRequest()

        then:
        config.getProject() >> 'my-project'
        config.getZone() >> ['my-zone']
        config.getRegion() >> ['my-region']
        config.getPreemptible() >> true
        // check staging script
        req.stagingScript == 'mkdir -p /work/dir; (alpha; beta; delta) 2>&1 > /work/dir/.command.log'
        req.diskSizeGb == null 
        // check main script
        req.mainScript == 'cd /work/dir; bash .command.run 2>&1 | tee -a .command.log'
        // check unstaging script
        req.unstagingScript.tokenize(';')[0] == 'cd /work/dir'
        req.unstagingScript.tokenize(';')[1] == ' [[ $GOOGLE_PIPELINE_FAILED == 1 || $NXF_DEBUG ]] && gsutil -m -q cp -R /google/ gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[2] == ' foo'
        req.unstagingScript.tokenize(';')[3] == ' bar'
        req.unstagingScript.tokenize(';')[4] == ' [[ -f .command.trace ]] && gsutil -m -q cp -R .command.trace gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[5] == ' gsutil -m -q cp -R .command.err gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[6] == ' gsutil -m -q cp -R .command.out gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[7] == ' gsutil -m -q cp -R .command.log gs://my-bucket/work/dir || true'
        req.unstagingScript.tokenize(';')[8] == ' gsutil -m -q cp -R .exitcode gs://my-bucket/work/dir || true'

        req.unstagingScript.tokenize(';').size() == 9
    }
    
    def 'should check if it is running'(){
        given:
        // -- task
        def task = Mock(TaskRun)
        task.name >> 'google-task'
        task.getWorkDir() >> Paths.get('/work/dir')

        // -- executor
        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)

        def operation = new Operation()

        // -- handler
        def handler = Spy(GooglePipelinesTaskHandler)
        handler.executor = executor
        handler.task = task
        handler.operation = operation

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
        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper:helper)

        // -- handler
        def handler = Spy(GooglePipelinesTaskHandler)
        handler.executor = executor
        handler.task = task

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
        def handler = [:] as GooglePipelinesTaskHandler
        expect:
        handler.getPipelineIdFromOp(operation) == '16737869387120678662'
    }

    def 'should get events from operation' () {
        given:
        def handler = [:] as GooglePipelinesTaskHandler

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

}
