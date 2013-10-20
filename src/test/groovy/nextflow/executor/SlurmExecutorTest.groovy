package nextflow.executor

import java.nio.file.Paths

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SlurmExecutorTest extends Specification {

    def testParseJob() {

        when:
        def exec = [:] as SlurmExecutor
        then:
        exec.parseJobId('Submitted batch job 10') == 10

    }

    def testKill() {

        given:
        def exec = [:] as SlurmExecutor
        expect:
        exec.killTaskCommand(123) == ['scancel','123']

    }

    def testGetSubmitCmdLine() {

        given:
        def base = Mock(BaseScript)
        def config = new TaskConfig(base)
        def script = Paths.get('/some/script.sh')
        def task = Mock(TaskRun)
        task.workDirectory >> Paths.get('/work/path')
        task.name >> 'task 555'
        def exec = [:] as SlurmExecutor
        exec.taskConfig = config

        when:
        config.maxDuration( Duration.create('1h') )
        config.clusterOptions = '-x -y -z'
        then:
        exec.getSubmitCommandLine(task,script).join(' ') == 'sbatch -D /work/path -J nf-task_555 -o /dev/null -t 01:00:00 -x -y -z script.sh'
    }

}
