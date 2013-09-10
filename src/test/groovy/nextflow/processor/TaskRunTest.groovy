package nextflow.processor

import spock.lang.Specification

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskRunTest extends Specification {

    def testGetInputsByType() {

        setup:
        def script = Mock(Script)
        def task = new TaskRun()
        task.setInput( new StdInParam(script), 'Hello' )
        task.setInput( new FileInParam(script, 'x'), 'file1' )
        task.setInput( new FileInParam(script, 'y'), 'file2' )
        task.setInput( new EnvInParam(script, 'z'), 'env' )


        when:
        def files = task.getInputsByType(FileInParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileInParam
        files.keySet()[1] instanceof FileInParam

        files.keySet()[0].name == 'x'
        files.keySet()[1].name == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }

    def testGetOutputsByType() {

        setup:
        def script = Mock(Script)
        def task = new TaskRun()
        task.setOutput( new FileOutParam(script, 'x'), 'file1' )
        task.setOutput( new FileOutParam(script, 'y'), 'file2' )
        task.setOutput( new StdOutParam(null), 'Hello' )


        when:
        def files = task.getOutputsByType(FileOutParam)
        then:
        files.size() == 2

        files.keySet()[0] instanceof FileOutParam
        files.keySet()[1] instanceof FileOutParam

        files.keySet()[0].name == 'x'
        files.keySet()[1].name == 'y'

        files.values()[0] == 'file1'
        files.values()[1] == 'file2'

    }



}