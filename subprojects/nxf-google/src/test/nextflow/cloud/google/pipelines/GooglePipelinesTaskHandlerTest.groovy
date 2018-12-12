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
import spock.lang.Specification

import java.nio.file.FileSystem
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.spi.FileSystemProvider

import com.google.api.services.genomics.v2alpha1.model.Operation
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CacheHelper

class GooglePipelinesTaskHandlerTest extends Specification {

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

        handler.operation == operation
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

    private Path mockPath(String path) {
        def provider = Mock(FileSystemProvider)
        provider.getScheme() >> 'gs'
        def fs = Mock(FileSystem)
        fs.provider() >> provider
        def uri = GroovyMock(URI)
        uri.toString() >> 'gs:/' + path

        def result = Mock(Path)
        result.toString() >> path
        result.toUriString() >> 'gs:/' + path
        result.getFileSystem() >> fs
        result.toUri() >> uri

        return result
    }

    def 'should create pipeline request' () {
        given:

        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)
        def config = Mock(GooglePipelinesConfiguration)
        def workDir = mockPath('/work/dir')

        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'

        def handler = new GooglePipelinesTaskHandler(
                pipelineConfiguration: config,
                executor: executor,
                task: task,
                instanceType: 'n1-1234'
        )

        when:
        def req = handler.createPipelineRequest()

        then:
        config.getProject() >> 'my-project'
        config.getZone() >> ['my-zone']
        config.getRegion() >> ['my-region']
        config.getPreemptible() >> true
        // chek request object
        req.instanceType == 'n1-1234'
        req.project == 'my-project'
        req.zone == ['my-zone']
        req.region == ['my-region']
        req.diskName == GooglePipelinesTaskHandler.diskName
        req.preemptible == true
        req.taskName == "nf-bad893071e9130b866d43a4fcabb95b6"
        req.containerImage == 'my/image'
        req.fileCopyImage == GooglePipelinesTaskHandler.fileCopyImage
        // check staging script
        req.stagingScript == 'mkdir -p /work/dir'
        // check main script
        req.mainScript == 'cd /work/dir; bash .command.run 2>&1 | tee -a .command.log'
        // check unstaging script
        req.unstagingScript.tokenize(';')[0] == '[[ $GOOGLE_PIPELINE_FAILED == 1 ]] && gsutil -m -q cp -c -P -r /google/ gs://work/dir || true'
        req.unstagingScript.tokenize(';')[1] == ' gsutil -m -q cp -c -P /work/dir/.command.err gs://work/dir || true'
        req.unstagingScript.tokenize(';')[2] == ' gsutil -m -q cp -c -P /work/dir/.command.out gs://work/dir || true'
        req.unstagingScript.tokenize(';')[3] == ' gsutil -m -q cp -c -P /work/dir/.command.log gs://work/dir || true'
        req.unstagingScript.tokenize(';')[4] == ' gsutil -m -q cp -c -P /work/dir/.exitcode gs://work/dir || true'
        req.unstagingScript.tokenize(';').size() == 5 
    }

    def 'should create pipeline request with stage and unstage commands' () {
        given:

        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper: helper)
        def config = Mock(GooglePipelinesConfiguration)
        def workDir = mockPath('/work/dir')

        def task = Mock(TaskRun)
        task.getName() >> 'foo'
        task.getWorkDir() >> workDir
        task.getHash() >> { CacheHelper.hasher('dummy').hash() }
        task.getContainer() >> 'my/image'

        def handler = new GooglePipelinesTaskHandler(
                pipelineConfiguration: config,
                executor: executor,
                instanceType: 'n1-1234',
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
        // check main script
        req.mainScript == 'cd /work/dir; bash .command.run 2>&1 | tee -a .command.log'
        // check unstaging script
        req.unstagingScript.tokenize(';')[0] == '[[ $GOOGLE_PIPELINE_FAILED == 1 ]] && gsutil -m -q cp -c -P -r /google/ gs://work/dir || true'
        req.unstagingScript.tokenize(';')[1] == ' foo'
        req.unstagingScript.tokenize(';')[2] == ' bar'
        req.unstagingScript.tokenize(';')[3] == ' gsutil -m -q cp -c -P /work/dir/.command.err gs://work/dir || true'
        req.unstagingScript.tokenize(';')[4] == ' gsutil -m -q cp -c -P /work/dir/.command.out gs://work/dir || true'
        req.unstagingScript.tokenize(';')[5] == ' gsutil -m -q cp -c -P /work/dir/.command.log gs://work/dir || true'
        req.unstagingScript.tokenize(';')[6] == ' gsutil -m -q cp -c -P /work/dir/.exitcode gs://work/dir || true'
        req.unstagingScript.tokenize(';').size() == 7
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


}
