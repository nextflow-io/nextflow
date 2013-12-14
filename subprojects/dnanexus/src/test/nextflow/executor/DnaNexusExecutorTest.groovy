package nextflow.executor
import java.nio.file.Files

import nextflow.fs.dx.api.DxApi
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DnaNexusExecutorTest extends Specification {

    def testCreateInputObject() {

        given:
        def exec = [:] as DnaNexusExecutor
        when:
        def obj  = exec.createInputObject( [a:1, b:2], 'dx_cc2.8xlarge' )
        then:
        obj.input == [a:1, b:2]
        obj.function == 'process'
        obj.systemRequirements.process.instanceType == 'dx_cc2.8xlarge'

    }


    def testHadlerSubmit() {

        given:
        def api = Mock(DxApi)
        def task = Mock(TaskRun)
        def exec = Mock(DnaNexusExecutor)
        def script = Mock(BaseScript)

        def params = [file1:'abc', file2: 'xxx']
        // note: the *instanceType* configured must used in the submit method
        def config = new TaskConfig(script)
        config.instanceType = 'dx_m1.super'
        // define the handler
        def handler = new DxTaskHandler(task, config, exec, params, api);

        // constraints
        and:
        1 * exec.createInputObject(params,'dx_m1.super') >> [function:'process']
        1 * api.jobNew(_) >> "job-xyz"

        when:
        handler.submit()
        then:
        handler.processJobId == "job-xyz"

    }


    def testKill() {

        given:
        def api = Mock(DxApi)
        def task = Mock(TaskRun)
        def exec = Mock(DnaNexusExecutor)
        def config = Mock(TaskConfig)
        def handler = new DxTaskHandler(task, config, exec, null, api);
        handler.processJobId = 'job-123'

        and:
        1 * api.jobTerminate('job-123')

        when:
        handler.kill()
        then:
        noExceptionThrown()

    }


    def testCheckStatus() {

        given:
        def api = Mock(DxApi)
        def task = Mock(TaskRun)
        def exec = Mock(DnaNexusExecutor)
        def config = Mock(TaskConfig)
        def handler = new DxTaskHandler(task, config, exec, null, api);
        handler.processJobId = 'job-312'

        and:
        1 * api.jobDescribe('job-312') >> [alpha:1, beta:2]

        when:
        handler.checkStatus() ==  [alpha:1, beta:2]
        then:
        noExceptionThrown()

    }


    def testHandleCheckIfRunning() {

        setup:
        def api = Mock(DxApi)
        def task = Mock(TaskRun)
        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)
        def handler = new DxTaskHandler(task, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'runnable'] }
        when:
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.SUBMITTED
        then:
        handler.checkIfRunning()

    }

    def testTaskHandlerCheckIfTerminated() {

        setup:
        def api = Mock(DxApi)
        def outFile = Files.createTempFile('testOutFile', null)
        def task = Mock(TaskRun)
        task.getCmdOutputFile() >> outFile
        def config = Mock(TaskConfig)
        def exec = Mock(DnaNexusExecutor)

        when:
        def handler = new DxTaskHandler(task, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'running'] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        !handler.checkIfCompleted()
        handler.status == TaskHandler.Status.RUNNING

        when:
        def task2 = new TaskRun()
        task2.workDirectory = Files.createTempDirectory('testHandler')
        task2.getCmdOutputFile().text = 'Task says Hola'
        handler = new DxTaskHandler(task2, config, exec, [:], api);
        handler.metaClass.checkStatus = { return [state:'done', output:[exit_code:33]] }
        handler.processJobId = '123'
        handler.status = TaskHandler.Status.RUNNING
        then:
        handler.checkIfCompleted()
        handler.status == TaskHandler.Status.COMPLETED
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