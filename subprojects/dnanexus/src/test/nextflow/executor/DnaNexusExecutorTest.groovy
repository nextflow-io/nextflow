package nextflow.executor
import java.nio.file.Files

import com.fasterxml.jackson.databind.JsonNode
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.util.DxHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {

    def testCreateInputObject() {

        when:
        JsonNode json = DnaNexusExecutor.createInputObject( [a:1, b:2], 'dx_cc2.8xlarge' )
        def obj = (Map)DxHelper.jsonToObj(json)
        then:
        obj.input == [a:1, b:2]
        obj.function == 'process'
        obj.systemRequirements.process.instanceType == 'dx_cc2.8xlarge'

    }

    def testInstanceTypeConfig( )  {

        setup:
        def script = Mock(BaseScript)
        def task = Mock(TaskRun)
        def config = new TaskConfig(script)
        config.instanceType = 'dx_cc2.8xlarge'

        when:
        def executor = new DnaNexusExecutor()
        executor.taskConfig = config
        def handler = executor.createTaskHandler(task)
        then:
        handler.taskConfig.instanceType == 'dx_cc2.8xlarge'


    }


    def testHandleCheckIfRunning() {

        setup:
        def task = Mock(TaskRun)
        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)
        def handler = new DxTaskHandler(task, config, exec, [:]);
        handler.metaClass.checkStatus = { return [state:'runnable'] }
        when:
        handler.processJobId = '123'
        then:
        handler.checkIfRunning()

    }

    def testTaskHandlerCheckIfTerminated() {

        setup:
        def outFile = Files.createTempFile('testOutFile', null)
        def task = Mock(TaskRun)
        task.getCmdOutputFile() >> outFile
        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)

        when:
        def handler = new DxTaskHandler(task, config, exec, [:]);
        handler.metaClass.checkStatus = { return [state:'running'] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        !handler.checkIfTerminated()
        handler.status == TaskHandler.Status.RUNNING

        when:
        def task2 = new TaskRun()
        task2.workDirectory = Files.createTempDirectory('testHandler')
        task2.getCmdOutputFile().text = 'Task says Hola'
        handler = new DxTaskHandler(task2, config, exec, [:]);
        handler.metaClass.checkStatus = { return [state:'done', output:[exit_code:33]] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        handler.checkIfTerminated()
        handler.status == TaskHandler.Status.TERMINATED
        task2.exitCode == 33
        task2.stdout == 'Task says Hola'


        cleanup:
        task2.workDirectory?.deleteDir()
    }



//    def void testFindFiles() {
//
//        setup:
//        def folder = Paths.get(URI.create('dxfs:///test_1/'))
//        Files.createDirectories( folder )
//
//        def target1 = folder.resolve('test_file1.txt')
//        def target2 = folder.resolve('test_file2.txt')
//        def target3 = folder.resolve('diff_name.fa')
//        Files.copy(new ByteArrayInputStream("Hello1".getBytes()), target1);
//        Files.copy(new ByteArrayInputStream("Hello2".getBytes()), target2);
//        Files.copy(new ByteArrayInputStream("Hello3".getBytes()), target3);
//
//        when:
//        Path folderUnderTest = Paths.get(URI.create('dxfs:///test_1/'))
//
//        List files = []
//        folderUnderTest.eachFileMatch(FileType.FILES, ~/test_.*/ ) { Path it -> files << it}
//
//        then:
//        files.size() == 2
//
//
//        cleanup:
//        folderUnderTest?.deleteDir()
//
//    }


}