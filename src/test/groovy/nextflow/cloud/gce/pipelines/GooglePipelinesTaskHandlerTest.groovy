package nextflow.cloud.gce.pipelines

import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Paths

import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus

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
        error.getMessage() == "No container specified for process $noContainerTaskRunner.name. Either specify the container to use in the process definition or with 'process.container' value in your config"
        !handler

    }

    def 'should construct correctly'() {
        when:
        def handler = new GooglePipelinesTaskHandler(stubTaskRunner ,stubExecutor,pipeConfig)

        then:
        handler.task.container == "testContainer"
        handler.taskName == "nf-task-$uuid-$handler.task.name"
        handler.taskInstanceName == "nf-task-$uuid-$handler.task.name-$handler.task.id"
    }

    def 'should submit a task'() {
        given:
        def handler = new GooglePipelinesTaskHandler(stubTaskRunner ,stubExecutor,pipeConfig)

        when:
        handler.submit()

        then:
        1 * stubExecutor.helper.createPipeline(_,_) >> new Pipeline()
        1* stubExecutor.helper.runPipeline(_,_) >> new Operation().setName("testOperation")
        handler.stagingCommands.size() == 3
        handler.unstagingCommands.size() == 4
    }

    def 'should check if it is running'(){
        given:
        // -- task
        def task = Mock(TaskRun)
        task.name >> 'gcp-task'
        task.getWorkDir() >> Paths.get('/work/dir')

        // -- executor
        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper)

        def OPERATION = new Operation()

        // -- handler
        def handler = Spy(GooglePipelinesTaskHandler)
        handler.executor = executor
        handler.task = task
        handler.operation = OPERATION

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

//        def stillRunning = handler.checkIfRunning()
//        def notRunning = handler.checkIfRunning()
//
//        then:
//        2 * stubExecutor.helper.checkOperationStatus(_) >>> [new Operation().setName("incomplete").setDone(false),new Operation().setName("complete").setDone(true)]
//
//        stillRunning
//        !notRunning
    }


    def 'should check if it is complete'() {
        given:
        // -- task
        def task = Mock(TaskRun)
        task.name >> 'gcp-task'
        task.getWorkDir() >> Paths.get('/work/dir')

        // -- executor
        def helper = Mock(GooglePipelinesHelper)
        def executor = new GooglePipelinesExecutor(helper)

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
