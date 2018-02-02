package nextflow.executor

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification

import java.nio.file.Paths

/**
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class TesTaskHandlerTest extends Specification {

    def 'should set resources' () {

        given:
        def executor = Mock(TesExecutor)
        def task = Mock(TaskRun)
        task.getName() >> 'tes-task'
        task.getWorkDir() >> Paths.get(".")
        task.getConfig() >> new TaskConfig(memory: '2GB', cpus: 4, disk: '10GB')
        def handler = new TesTaskHandler(task, executor)


        when:
        def t = handler.newTesTask()

        then:
        t.getResources().cpuCores == 4
        t.getResources().ramGb == 2
        t.getResources().diskGb == 10
        t.getResources().preemptible == null
        t.getResources().zones == null

    }
}
