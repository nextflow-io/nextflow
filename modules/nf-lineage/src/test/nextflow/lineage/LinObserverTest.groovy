/*
 * Copyright 2013-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.lineage

import nextflow.extension.FilesEx
import nextflow.lineage.exception.OutputRelativePathException

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.file.FileHolder
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.TaskOutput
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.config.LineageConfig
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import nextflow.script.TokenVar
import nextflow.script.WorkflowMetadata
import nextflow.script.params.EnvOutParam
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.script.params.StdInParam
import nextflow.script.params.StdOutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.trace.event.WorkflowOutputEvent
import nextflow.util.CacheHelper
import nextflow.util.PathNormalizer
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

import static nextflow.lineage.fs.LinPath.*

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LinObserverTest extends Specification {
    @Shared
    Path lidFolder = Files.createTempDirectory("lid")
    def cleanupSpec(){
            lidFolder.deleteDir()
    }

    def 'should normalize paths' (){
        given:
        def folder = Files.createTempDirectory('test')
        def workDir = folder.resolve("workDir")
        def projectDir = folder.resolve("projectDir")
        def metadata = Mock(WorkflowMetadata){
            getRepository() >> "https://nextflow.io/nf-test/"
            getCommitId() >> "123456"
            getScriptId() >> "78910"
            getProjectDir() >> projectDir
            getWorkDir() >> workDir
        }
        def params = [path: workDir.resolve("path/file.txt"), sequence: projectDir.resolve("file2.txt").toString(), value: 12]
        when:
        def results = LinObserver.getNormalizedParams(params, new PathNormalizer(metadata))
        then:
        results.size() == 3
        results.get(0).name == "path"
        results.get(0).type == Path.simpleName
        results.get(0).value == "work/path/file.txt"
        results.get(1).name == "sequence"
        results.get(1).type == "String"
        results.get(1).value == projectDir.resolve("file2.txt").toString()
        results.get(2).name == "value"
        results.get(2).type == "Integer"
        results.get(2).value == 12

        cleanup:
        ScriptMeta.reset()
        folder?.deleteDir()
    }
    def 'should collect script files' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        def config = [workflow:[lineage:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultLinStore();
        def uniqueId = UUID.randomUUID()
        def scriptFile = folder.resolve("main.nf")
        def module1 = folder.resolve("script1.nf"); module1.text = 'hola'
        def module2 = folder.resolve("script2.nf"); module2.text = 'world'
        and:

        def metadata = Mock(WorkflowMetadata){
            getRepository() >> "https://nextflow.io/nf-test/"
            getCommitId() >> "123456"
            getScriptId() >> "78910"
            getScriptFile() >> scriptFile
            getProjectDir() >> folder.resolve("projectDir")
            getWorkDir() >> folder.resolve("workDir")
        }
        def session = Mock(Session) {
            getConfig() >> config
            getUniqueId() >> uniqueId
            getRunName() >> "test_run"
            getWorkflowMetadata() >> metadata
            getParams() >> new ScriptBinding.ParamsMap()
        }
        store.open(LineageConfig.create(session))
        def observer = Spy(new LinObserver(session, store))

        when:
        def files = observer.collectScriptDataPaths(new PathNormalizer(metadata))
        then:
        observer.allScriptFiles() >> [ scriptFile, module1, module2 ]
        and:
        files.size() == 3
        and:
        files[0].path == "file://${scriptFile.toString()}"
        files[0].checksum == new Checksum("78910", "nextflow", "standard")
        and:
        files[1].path == "file://$module1"
        files[1].checksum == Checksum.ofNextflow(module1.text)
        and:
        files[2].path == "file://$module2"
        files[2].checksum == Checksum.ofNextflow(module2.text)

        cleanup:
        ScriptMeta.reset()
        folder?.deleteDir()
    }

    def 'should save workflow' (){
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage:[enabled: true, store:[location:folder.toString()]]]
        def store = new DefaultLinStore();
        def uniqueId = UUID.randomUUID()
        def scriptFile = folder.resolve("main.nf")
        def map = [
            repository: "https://nextflow.io/nf-test/",
            commitId: "123456",
            scriptId: "78910",
            scriptFile: scriptFile,
            projectDir: folder.resolve("projectDir"),
            revision: "main",
            projectName: "nextflow.io/nf-test",
            workDir: folder.resolve("workDir")
        ]
        def metadata = Mock(WorkflowMetadata){
            getRepository() >> map.repository
            getCommitId() >> map.commitId
            getScriptId() >> map.scriptId
            getScriptFile() >> map.scriptFile
            getProjectDir() >> map.projectDir
            getRevision() >> map.revision
            getProjectName() >> map.projectName
            getWorkDir() >> map.workDir
            toMap() >> map
        }
        def session = Mock(Session) {
            getConfig() >> config
            getUniqueId() >> uniqueId
            getRunName() >> "test_run"
            getWorkflowMetadata() >> metadata
            getParams() >> new ScriptBinding.ParamsMap()
        }
        store.open(LineageConfig.create(session))
        def observer = new LinObserver(session, store)
        def mainScript = new DataPath("file://${scriptFile.toString()}", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], map.repository, map.commitId)
        def workflowRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [], config, map)
        when:
        observer.onFlowCreate(session)
        observer.onFlowBegin()
        then:
        folder.resolve("${observer.executionHash}/.data.json").text == new LinEncoder().encode(workflowRun)

        cleanup:
        folder?.deleteDir()
    }

    @Unroll
    def 'should get parameter type' () {
        expect:
        LinObserver.getParameterType(PARAM) == STRING
        where:
        PARAM                                           | STRING
        null                                            | null
        new FileInParam(null, [])                       | "path"
        new ValueOutParam(null, [])                     | "val"
        new EnvOutParam(null, [])                       | "env"
        new StdInParam(null, [])                        | "stdin"
        new StdOutParam(null, [])                       | "stdout"
        Path.of("test")                                 | "Path"
        ["test"]                                        | "Collection"
        [key:"value"]                                   | "Map"
    }

    def 'should save task run' () {
        given:
        def folder = Files.createTempDirectory('test').toRealPath()
        def config = [workflow:[lineage:[enabled: true, store:[location:folder.toString()]]]]
        def uniqueId = UUID.randomUUID()
        def workDir = folder.resolve("work")
        def session = Mock(Session) {
            getConfig()>>config
            getUniqueId()>>uniqueId
            getRunName()>>"test_run"
            getWorkDir() >> workDir
        }
        def metadata = Mock(WorkflowMetadata){
            getRepository() >> "https://nextflow.io/nf-test/"
            getCommitId() >> "123456"
            getScriptId() >> "78910"
            getProjectDir() >> folder.resolve("projectDir")
            getWorkDir() >> workDir
        }
        and:
        def store = new DefaultLinStore();
        store.open(LineageConfig.create(session))
        and:
        def observer = Spy(new LinObserver(session, store))
        def normalizer = new PathNormalizer(metadata)
        observer.executionHash = "hash"
        observer.normalizer = normalizer
        and:
        def hash = HashCode.fromString("1234567890")
        def taskWd = workDir.resolve('12/34567890')
        Files.createDirectories(taskWd)
        and:
        observer.getTaskGlobalVars(_) >> [:]
        observer.getTaskBinEntries(_) >> []

        and: 'Task Inputs'
        def inputs = new LinkedHashMap<InParam, Object>()
        // File from task
        inputs.put(new FileInParam(null, []).bind("file1"), [new FileHolder(workDir.resolve('78/567890/file1.txt'))])
        // Normal file
        def file = folder.resolve("file2.txt")
        file.text = "this is a test file"
        def fileHash = CacheHelper.hasher(file).hash().toString()
        inputs.put(new FileInParam(null, []).bind("file2"), [new FileHolder(file)])
        //Value input
        inputs.put(new ValueInParam(null, []).bind("id"), "value")

        and: 'Task Outputs'
        def outputs = new LinkedHashMap<OutParam, Object>()
        // Single Path output
        def outFile1 = taskWd.resolve('fileOut1.txt')
        outFile1.text = 'some data'
        def fileHash1 = CacheHelper.hasher(outFile1).hash().toString()
        def attrs1 = Files.readAttributes(outFile1, BasicFileAttributes)
        outputs.put(new FileOutParam(null, []).bind(new TokenVar("file1")), outFile1)
        // Collection Path output
        def outFile2 = taskWd.resolve('fileOut2.txt')
        outFile2.text = 'some other data'
        def fileHash2 = CacheHelper.hasher(outFile2).hash().toString()
        def attrs2 = Files.readAttributes(outFile2, BasicFileAttributes)
        outputs.put(new FileOutParam(null, []).bind(new TokenVar("file2")), [outFile2])
        outputs.put(new ValueOutParam(null, []).bind(new TokenVar("id")), "value")

        and: 'Task description'
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getSource() >> 'echo task source'
            getScript() >> 'this is the script'
            getInputs() >> inputs
            getOutputs() >> outputs
            getWorkDir() >> taskWd
        }
        def handler = Mock(TaskHandler){
            getTask() >> task
        }

        and: 'Expected LID objects'
        def sourceHash = CacheHelper.hasher('echo task source').hash().toString()
        def script = 'this is the script'
        def taskDescription = new nextflow.lineage.model.v1beta1.TaskRun(uniqueId.toString(), "foo",
            new Checksum(sourceHash, "nextflow", "standard"),
            script,
            [
                new Parameter("path", "file1", ['lid://78567890/file1.txt']),
                new Parameter("path", "file2", [[path: normalizer.normalizePath(file), checksum: [value:fileHash, algorithm: "nextflow", mode:  "standard"]]]),
                new Parameter("val", "id", "value")
            ], null, null, null, null, [:], [], "lid://hash")
        def dataOutput1 = new FileOutput(outFile1.toString(), new Checksum(fileHash1, "nextflow", "standard"),
            "lid://1234567890", "lid://hash", "lid://1234567890", attrs1.size(), LinUtils.toDate(attrs1.creationTime()), LinUtils.toDate(attrs1.lastModifiedTime()) )
        def dataOutput2 = new FileOutput(outFile2.toString(), new Checksum(fileHash2, "nextflow", "standard"),
            "lid://1234567890", "lid://hash", "lid://1234567890", attrs2.size(), LinUtils.toDate(attrs2.creationTime()), LinUtils.toDate(attrs2.lastModifiedTime()) )

        when:
        observer.onTaskComplete(new TaskEvent(handler, null))
        def taskRunResult = store.load("${hash.toString()}")
        def dataOutputResult1 = store.load("${hash}/fileOut1.txt") as FileOutput
        def dataOutputResult2 = store.load("${hash}/fileOut2.txt") as FileOutput
        def taskOutputsResult = store.load("${hash}#output") as TaskOutput
        then:
        taskRunResult == taskDescription
        dataOutputResult1 == dataOutput1
        dataOutputResult2 == dataOutput2
        taskOutputsResult.taskRun == "lid://1234567890"
        taskOutputsResult.workflowRun == "lid://hash"
        taskOutputsResult.output.size() == 3
        taskOutputsResult.output.get(0).type == "path"
        taskOutputsResult.output.get(0).name == "file1"
        taskOutputsResult.output.get(0).value == "lid://1234567890/fileOut1.txt"
        taskOutputsResult.output.get(1).type == "path"
        taskOutputsResult.output.get(1).name == "file2"
        taskOutputsResult.output.get(1).value == ["lid://1234567890/fileOut2.txt"]
        taskOutputsResult.output.get(2).type == "val"
        taskOutputsResult.output.get(2).name == "id"
        taskOutputsResult.output.get(2).value == "value"

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task data output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage:[enabled: true, store:[location:folder.toString()]]]
        def store = new DefaultLinStore();
        def session = Mock(Session) {
            getConfig()>>config
        }
        store.open(LineageConfig.create(session))
        def observer = Spy(new LinObserver(session, store))
        observer.executionHash = "hash"
        and:
        def workDir = folder.resolve('12/34567890')
        Files.createDirectories(workDir)
        and:
        def outFile = workDir.resolve('foo/bar/file.bam')
        Files.createDirectories(outFile.parent)
        outFile.text = 'some data'
        def fileHash = CacheHelper.hasher(outFile).hash().toString()
        and:
        def hash = HashCode.fromInt(123456789)
        and:
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getWorkDir() >> workDir
        }
        and:
        def attrs = Files.readAttributes(outFile, BasicFileAttributes)
        def output = new FileOutput(outFile.toString(), new Checksum(fileHash, "nextflow", "standard"),
            "lid://15cd5b07", "lid://hash", "lid://15cd5b07", attrs.size(), LinUtils.toDate(attrs.creationTime()), LinUtils.toDate(attrs.lastModifiedTime()) )
        and:
        observer.readAttributes(outFile) >> attrs

        when:
        observer.storeTaskOutput(task, outFile)
        then:
        folder.resolve("${hash}/foo/bar/file.bam/.data.json").text == new LinEncoder().encode(output)

        cleanup:
        folder?.deleteDir()
    }

    def 'should relativise task output dirs' (){
        when:
        def config = [workflow:[lineage:[enabled: true, store:[location:lidFolder.toString()]]]]
        def store = new DefaultLinStore();
        def session = Mock(Session) {
            getConfig()>>config
        }
        def hash = HashCode.fromInt(123456789)
        def taskConfig = Mock(TaskConfig){
            getStoreDir() >> STORE_DIR
        }
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getWorkDir() >> WORK_DIR
            getConfig() >> taskConfig
        }
        store.open(LineageConfig.create(session))
        def observer = new LinObserver(session, store)
        then:
        observer.getTaskRelative(task, PATH) == EXPECTED
        where:
        WORK_DIR                            | STORE_DIR                         | PATH                                          | EXPECTED
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('/path/to/work/12/3456789/relative')  | "relative"
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('/path/to/storeDir/relative')         | "relative"
        Path.of('work/12/3456789')          | Path.of('storeDir')               | Path.of('work/12/3456789/relative')           | "relative"
        Path.of('work/12/3456789')          | Path.of('storeDir')               | Path.of('storeDir/relative')                  | "relative"
        Path.of('work/12/3456789')          | Path.of('storeDir')               | Path.of('results/relative')                   | "results/relative"
        Path.of('/path/to/work/12/3456789') | Path.of('storeDir')               | Path.of('./relative')                         | "relative"
    }

    @Unroll
    def 'should return exception when relativize task output dirs'() {
        when:
            def config = [workflow:[lineage:[enabled: true, store:[location:lidFolder.toString()]]]]
            def store = new DefaultLinStore();
            def session = Mock(Session) {
                getConfig()>>config
            }
        def hash = HashCode.fromInt(123456789)
        def taskConfig = Mock(TaskConfig){
            getStoreDir() >> STORE_DIR
        }
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getWorkDir() >> WORK_DIR
            getConfig() >> taskConfig
        }
        store.open(LineageConfig.create(session))
        def observer = new LinObserver(session, store)
        observer.getTaskRelative(task, PATH)
        then:
            def e = thrown(IllegalArgumentException)
            e.message == "Cannot access the relative path for output '$PATH' and task '${task.name}'".toString()

        where:
        WORK_DIR                            | STORE_DIR                         | PATH
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('/another/path/relative')
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('../path/to/storeDir/relative')
    }

    def 'should relativize workflow output dirs' (){
        when:
            def config = [workflow:[lineage:[enabled: true, store:[location:lidFolder.toString()]]]]
            def store = new DefaultLinStore();
            def session = Mock(Session) {
                getOutputDir()>>OUTPUT_DIR
                getConfig()>>config
            }
            store.open(LineageConfig.create(session))
            def observer = new LinObserver(session, store)
        then:
            observer.getWorkflowRelative(PATH) == EXPECTED
        where:
        OUTPUT_DIR                      | PATH                                  | EXPECTED
        Path.of('/path/to/outDir')      | Path.of('/path/to/outDir/relative')   | "relative"
        Path.of('outDir')               | Path.of('outDir/relative')            | "relative"
        Path.of('/path/to/outDir')      | Path.of('results/relative')           | "results/relative"
        Path.of('/path/to/outDir')      | Path.of('./relative')                 | "relative"
    }

    @Unroll
    def 'should return exception when relativize workflow output dirs' (){
        when:
            def config = [workflow:[lineage:[enabled: true, store:[location:lidFolder.toString()]]]]
            def store = new DefaultLinStore();
            def session = Mock(Session) {
                getOutputDir()>>OUTPUT_DIR
                getConfig()>>config
            }
            def observer = new LinObserver(session, store)
            observer.getWorkflowRelative(PATH)
        then:
            thrown(OutputRelativePathException)
        where:
        OUTPUT_DIR                      | PATH                                  | EXPECTED
        Path.of('/path/to/outDir')      | Path.of('/another/path/')             | "relative"
        Path.of('/path/to/outDir')      | Path.of('../relative')                | "relative"
    }

    def 'should save workflow output'() {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage:[enabled: true, store:[location:folder.toString()]]]
        def store = new DefaultLinStore();
        def outputDir = folder.resolve('results')
        def uniqueId = UUID.randomUUID()
        def scriptFile = folder.resolve("main.nf")
        def workDir= folder.resolve("work")
        def metadata = Mock(WorkflowMetadata){
            getRepository() >> "https://nextflow.io/nf-test/"
            getCommitId() >> "123456"
            getScriptId() >> "78910"
            getScriptFile() >> scriptFile
            getProjectDir() >> folder.resolve("projectDir")
            getWorkDir() >> workDir
        }
        def session = Mock(Session) {
            getConfig()>>config
            getOutputDir()>>outputDir
            getWorkDir() >> workDir
            getWorkflowMetadata()>>metadata
            getUniqueId()>>uniqueId
            getRunName()>>"test_run"
            getParams() >> new ScriptBinding.ParamsMap()
        }
        store.open(LineageConfig.create(session))
        def observer = new LinObserver(session, store)
        def encoder = new LinEncoder()

        when: 'Starting workflow'
            observer.onFlowCreate(session)
            observer.onFlowBegin()
        then: 'History file should contain execution hash'
            def lid = LinHistoryRecord.parse(folder.resolve(".history/${observer.executionHash}").text)
            lid.runLid == asUriString(observer.executionHash)
            lid.sessionId == uniqueId
            lid.runName == "test_run"

        when: ' publish output with source file'
            def outFile1 = outputDir.resolve('foo/file.bam')
            Files.createDirectories(outFile1.parent)
            outFile1.text = 'some data1'
            def sourceFile1 = workDir.resolve('12/3987/file.bam')
            Files.createDirectories(sourceFile1.parent)
            sourceFile1.text = 'some data1'
            observer.onFilePublish(new FilePublishEvent(sourceFile1, outFile1))
            observer.onWorkflowOutput(new WorkflowOutputEvent("a", outFile1))

        then: 'check file 1 output metadata in lid store'
            def attrs1 = Files.readAttributes(outFile1, BasicFileAttributes)
            def fileHash1 = CacheHelper.hasher(outFile1).hash().toString()
            def output1 = new FileOutput(outFile1.toString(), new Checksum(fileHash1, "nextflow", "standard"),
                "lid://123987/file.bam", "$LID_PROT${observer.executionHash}", null,
                attrs1.size(), LinUtils.toDate(attrs1.creationTime()), LinUtils.toDate(attrs1.lastModifiedTime()) )
            folder.resolve("${observer.executionHash}/foo/file.bam/.data.json").text == encoder.encode(output1)

        when: 'publish without source path'
        def outFile2 = outputDir.resolve('foo/file2.bam')
            Files.createDirectories(outFile2.parent)
            outFile2.text = 'some data2'
            def attrs2 = Files.readAttributes(outFile2, BasicFileAttributes)
            def fileHash2 = CacheHelper.hasher(outFile2).hash().toString()
            observer.onFilePublish(new FilePublishEvent(null, outFile2))
            observer.onWorkflowOutput(new WorkflowOutputEvent("b", outFile2))
        then: 'Check outFile2 metadata in lid store'
            def output2 = new FileOutput(outFile2.toString(), new Checksum(fileHash2, "nextflow", "standard"),
                "lid://${observer.executionHash}" , "lid://${observer.executionHash}", null,
                attrs2.size(), LinUtils.toDate(attrs2.creationTime()), LinUtils.toDate(attrs2.lastModifiedTime()) )
            folder.resolve("${observer.executionHash}/foo/file2.bam/.data.json").text == encoder.encode(output2)

        when: 'Workflow complete'
            observer.onFlowComplete()
        then: 'Check WorkflowOutput is written in the lid store'
            def resultsRetrieved = store.load("${observer.executionHash}#output") as WorkflowOutput
            resultsRetrieved.output == [new Parameter(Path.simpleName, "a", "lid://${observer.executionHash}/foo/file.bam"), new Parameter(Path.simpleName, "b", "lid://${observer.executionHash}/foo/file2.bam")]

        cleanup:
            folder?.deleteDir()
    }

    def 'should not save workflow output entry when no outputs'() {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage: [enabled: true, store: [location: folder.toString()]]]
        def store = new DefaultLinStore();
        def outputDir = folder.resolve('results')
        def uniqueId = UUID.randomUUID()
        def scriptFile = folder.resolve("main.nf")
        def workDir = folder.resolve("work")
        def metadata = Mock(WorkflowMetadata) {
            getRepository() >> "https://nextflow.io/nf-test/"
            getCommitId() >> "123456"
            getScriptId() >> "78910"
            getScriptFile() >> scriptFile
            getProjectDir() >> folder.resolve("projectDir")
            getWorkDir() >> workDir
        }
        def session = Mock(Session) {
            getConfig() >> config
            getOutputDir() >> outputDir
            getWorkDir() >> workDir
            getWorkflowMetadata() >> metadata
            getUniqueId() >> uniqueId
            getRunName() >> "test_run"
            getParams() >> new ScriptBinding.ParamsMap()
        }
        store.open(LineageConfig.create(session))
        def observer = new LinObserver(session, store)

        when:
        observer.onFlowCreate(session)
        observer.onFlowBegin()
        observer.onFlowComplete()
        def resultFile = folder.resolve("${observer.executionHash}#output")
        then:
        !resultFile.exists()
    }
}
