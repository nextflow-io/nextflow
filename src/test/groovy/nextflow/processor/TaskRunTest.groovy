package nextflow.processor

import java.nio.file.Paths

import nextflow.script.EnvInParam
import nextflow.script.FileInParam
import nextflow.script.FileOutParam
import nextflow.script.FileSharedParam
import nextflow.script.StdInParam
import nextflow.script.StdOutParam
import nextflow.script.TokenVar
import nextflow.script.ValueInParam
import nextflow.script.ValueOutParam
import spock.lang.Specification

/**
 *
 *  @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

class TaskRunTest extends Specification {

    def testGetInputsByType() {

        setup:
        def binding = new Binding('x': 1, 'y': 2)
        def task = new TaskRun()
        def list = []

        task.setInput( new StdInParam(binding,list) )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('x')), 'file1' )
        task.setInput( new FileInParam(binding, list).bind(new TokenVar('y')), 'file2' )
        task.setInput( new EnvInParam(binding, list).bind('z'), 'env' )


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
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        task.setOutput( new FileOutParam(binding, list).bind('x'), 'file1' )
        task.setOutput( new FileOutParam(binding, list).bind('y'), 'file2' )
        task.setOutput( new StdOutParam(binding, list), 'Hello' )


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
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        def x = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def y = new FileInParam(binding, list).bind('y')
        def z = new FileSharedParam(binding, list).bind('z')

        task.setInput(x, 1)
        task.setInput(y, [ new FileHolder(Paths.get('file_y_1')) ])
        task.setInput(z, [ new FileHolder(Paths.get('file_z_1')), new FileHolder(Paths.get('file_z_2')) ])

        expect:
        task.getInputFiles().size() == 2
        task.stagedInputs.size() == 3

    }


    def testGetOutputFilesNames() {

        setup:
        def binding = new Binding()
        def task = new TaskRun()
        def list = []

        when:
        def i1 = new ValueInParam(binding, list).bind( new TokenVar('x') )
        def s1 = new FileSharedParam(binding, list).bind( new TokenVar('y') )
        def o1 = new FileOutParam(binding,list).bind('file_out.alpha')
        def o2 = new ValueOutParam(binding,list).bind( 'x' )
        def o3 = new FileOutParam(binding,list).bind('file_out.beta')

        task.setInput(i1, 'Hello' )
        task.setInput(s1, [ new FileHolder(Paths.get('file_shared.delta')) ])
        task.setOutput(o1)
        task.setOutput(o2)
        task.setOutput(o3)

        then:
        task.getOutputFilesNames() == ['file_out.alpha', 'file_out.beta', 'file_shared.delta']
    }



    def testDumpStdout() {

        setup:
        def task = new TaskRun()
        task.stdout = """
            1
            2
            3
            4
            5
            6
            7
            8
            9
            """.stripIndent()

        when:
        def lines = []
        def count = task.dumpStdout(lines, 5)

        then:
        count == 5
        lines.size() == 5
        lines[0] == '  '
        lines[1] == '  1'
        lines[2] == '  2'
        lines[3] == '  3'
        lines[4] == '  4'



    }


}