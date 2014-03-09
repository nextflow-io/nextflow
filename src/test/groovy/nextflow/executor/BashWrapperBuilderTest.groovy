package nextflow.executor

import nextflow.processor.FileHolder
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.FileInParam
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BashWrapperBuilderTest extends Specification {

    def 'test changeToScratchDir' () {

        setup:
        def builder = [:] as BashWrapperBuilder

        expect:
        builder.changeToScratchDirectory() == null

        when:
        builder.scratch = true
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NF_SCRATCH'

        when:
        builder.scratch = '$SOME_DIR'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$SOME_DIR && cd $NF_SCRATCH'

        when:
        builder.scratch = '$SOME_DIR'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$SOME_DIR && cd $NF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

        when:
        builder.scratch = '/my/temp'
        then:
        builder.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

    }


    def testDocketMounts() {

        when:
        def task = new TaskRun()
        def builder = new BashWrapperBuilder(task)
        then:
        builder.getDockerMounts() == '-v $PWD:$PWD'


        when:
        task = new TaskRun()
        task.inputs[ new FileInParam(null,[]) ] = [FileHolder.get('/home/data/sequences', 'file.txt')]
        task.inputs[ new FileInParam(null,[]) ] = [FileHolder.get('/home/data/file1','seq_1.fa'), FileHolder.get('/home/data/file2','seq_2.fa'), FileHolder.get('/home/data/file3','seq_3.fa') ]
        builder = new BashWrapperBuilder(task)
        then:
        builder.getDockerMounts() == '-v /home/data:/home/data -v $PWD:$PWD'

    }

    def testBashWrapperConstructor() {

        setup:
        def process = Mock(TaskProcessor)
        process.getProcessEnvironment() >> [VAR_1:'X', VAR_2: 'Y']

        def task = new TaskRun()
        task.stdin = 'Hello'
        task.processor = process

        def config = new TaskConfig(Mock(BaseScript), [dockerize:'fedora', scratch:true])

        when:
        def bash = new BashWrapperBuilder(task, config)
        then:
        bash.input == 'Hello'
        bash.scratch == true
        bash.dockerContainerName == 'fedora'
        bash.environment == [VAR_1:'X', VAR_2: 'Y']

    }

}
