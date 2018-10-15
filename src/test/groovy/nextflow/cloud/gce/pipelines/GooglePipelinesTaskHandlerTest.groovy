package nextflow.cloud.gce.pipelines


import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Paths

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
        getHelper() >> new GooglePipelinesHelper()
    }


    TaskRun stubTaskRunner = GroovyStub {
        getName() >> "testName"
        getId() >> new TaskId(12345)
        getContainer() >> "testContainer"
        getScript() >> "echo testScript"
        getWorkDir() >> Paths.get("/")
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
}