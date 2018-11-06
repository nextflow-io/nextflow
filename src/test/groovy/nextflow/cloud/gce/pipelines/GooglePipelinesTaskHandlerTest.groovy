package nextflow.cloud.gce.pipelines

import com.google.api.services.genomics.v2alpha1.model.Operation
import com.google.api.services.genomics.v2alpha1.model.Pipeline
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import spock.lang.Shared
import spock.lang.Specification

class GooglePipelinesTaskHandlerTest extends Specification {

    @Shared
    GooglePipelinesConfiguration pipeConfig = new GooglePipelinesConfiguration("testProject","testZone","instanceType")

    @Shared
    UUID uuid = new UUID(4,4)

    @Shared
    Session stubSession = Stub {
        getUniqueId() >> uuid
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
        error.getMessage() == "No container is specified for process $noContainerTaskRunner.name . Either specify the container to use in the process definition or with 'process.container' value in your config"
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
        def handler = new GooglePipelinesTaskHandler(stubTaskRunner,stubExecutor,pipeConfig)

        when:
        def stillRunning = handler.checkIfRunning()
        def notRunning = handler.checkIfRunning()

        then:
        2 * stubExecutor.helper.checkOperationStatus(_) >>> [new Operation().setName("incomplete").setDone(false),new Operation().setName("complete").setDone(true)]

        stillRunning
        !notRunning
    }

    def 'should check if it is complete'() {
        given:
        def handler = new GooglePipelinesTaskHandler(stubTaskRunner,stubExecutor,pipeConfig)

        when:
        def notComplete = handler.checkIfCompleted()
        def complete = handler.checkIfCompleted()

        then:
        2 * stubExecutor.helper.checkOperationStatus(_) >>> [new Operation().setName("incomplete").setDone(false),new Operation().setName("complete").setDone(true)]

        !notComplete
        complete
    }


}
