package nextflow.processor

import java.nio.file.Paths

import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.ScriptVar
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.ValueInParam
import spock.lang.Specification

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskRunTest extends Specification {

    def testGetInputsByType() {

        setup:
        def script = Mock(Script)
        script.getBinding() >> { new Binding('x': 1, 'y': 2) }
        def task = new TaskRun()
        task.setInput( new StdInParam(script,'Hello') )
        task.setInput( new FileInParam(script, new ScriptVar('x')), 'file1' )
        task.setInput( new FileInParam(script, new ScriptVar('y')), 'file2' )
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

    def testGetInputFiles() {

        setup:
        def script = Mock(Script)
        def task = new TaskRun()

        def x = new ValueInParam(script, 'x')
        def y = new FileInParam(script, 'y')
        def z = new FileSharedParam(script, 'z')

        task.setInput(x, 1)
        task.setInput(y, [ new FileHolder(Paths.get('file_y_1')) ])
        task.setInput(z, [ new FileHolder(Paths.get('file_z_1')), new FileHolder(Paths.get('file_z_2')) ])

        expect:
        task.getInputFiles().size() == 2
        task.stagedInputs.size() == 3

    }



}