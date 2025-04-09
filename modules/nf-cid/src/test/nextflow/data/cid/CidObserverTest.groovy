/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.data.cid

import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.TaskOutputs
import nextflow.file.FileHolder
import nextflow.processor.TaskHandler
import nextflow.script.TokenVar
import nextflow.script.params.FileInParam
import nextflow.script.params.FileOutParam
import nextflow.script.params.InParam
import nextflow.script.params.OutParam
import nextflow.script.params.ValueInParam
import nextflow.script.params.ValueOutParam

import static nextflow.data.cid.fs.CidPath.*

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowOutputs
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.cid.serde.CidEncoder
import nextflow.data.config.DataConfig
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptMeta
import nextflow.script.WorkflowMetadata
import nextflow.util.CacheHelper
import nextflow.util.PathNormalizer
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidObserverTest extends Specification {
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
        def results = CidObserver.getNormalizedParams(params, new PathNormalizer(metadata))
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
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultCidStore();
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
        store.open(DataConfig.create(session))
        def observer = Spy(new CidObserver(session, store))

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
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultCidStore();
        def uniqueId = UUID.randomUUID()
        def scriptFile = folder.resolve("main.nf")
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
        store.open(DataConfig.create(session))
        def observer = new CidObserver(session, store)
        def mainScript = new DataPath("file://${scriptFile.toString()}", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript],"https://nextflow.io/nf-test/", "123456" )
        def workflowRun = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [])
        when:
        observer.onFlowCreate(session)
        observer.onFlowBegin()
        then:
        folder.resolve(".meta/${observer.executionHash}/.data.json").text == new CidEncoder().encode(workflowRun)

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task run' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
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
        def store = new DefaultCidStore();
        store.open(DataConfig.create(session))
        and:
        def observer = new CidObserver(session, store)
        def normalizer = new PathNormalizer(metadata)
        observer.executionHash = "hash"
        observer.normalizer = normalizer
        and:
        def hash = HashCode.fromString("1234567890")
        def taskWd = workDir.resolve('12/34567890')
        Files.createDirectories(taskWd)
        and:
        def processor = Mock(TaskProcessor){
            getTaskGlobalVars(_) >> [:]
            getTaskBinEntries(_) >> []
        }

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
            getProcessor() >> processor
            getSource() >> 'echo task source'
            getScript() >> 'this is the script'
            getInputs() >> inputs
            getOutputs() >> outputs
            getWorkDir() >> taskWd
        }
        def handler = Mock(TaskHandler){
            getTask() >> task
        }

        and: 'Expected CID objects'
        def sourceHash = CacheHelper.hasher('echo task source').hash().toString()
        def scriptHash = CacheHelper.hasher('this is the script').hash().toString()
        def taskDescription = new nextflow.data.cid.model.TaskRun(uniqueId.toString(), "foo",
            new Checksum(sourceHash, "nextflow", "standard"),
            new Checksum(scriptHash, "nextflow", "standard"),
            [
                new Parameter(FileInParam.simpleName, "file1", ['cid://78567890/outputs/file1.txt']),
                new Parameter(FileInParam.simpleName, "file2", [[path: normalizer.normalizePath(file), checksum: [value:fileHash, algorithm: "nextflow", mode:  "standard"]]]),
                new Parameter(ValueInParam.simpleName, "id", "value")
            ], null, null, null, null, [:], [], "cid://hash", null)
        def dataOutput1 = new DataOutput(outFile1.toString(), new Checksum(fileHash1, "nextflow", "standard"),
            "cid://1234567890", "cid://1234567890", attrs1.size(), CidUtils.toDate(attrs1.creationTime()), CidUtils.toDate(attrs1.lastModifiedTime()) )
        def dataOutput2 = new DataOutput(outFile2.toString(), new Checksum(fileHash2, "nextflow", "standard"),
            "cid://1234567890", "cid://1234567890", attrs2.size(), CidUtils.toDate(attrs2.creationTime()), CidUtils.toDate(attrs2.lastModifiedTime()) )

        when:
        observer.onProcessComplete(handler, null )
        def taskRunResult = store.load("${hash.toString()}")
        def dataOutputResult1 = store.load("${hash}/outputs/fileOut1.txt") as DataOutput
        def dataOutputResult2 = store.load("${hash}/outputs/fileOut2.txt") as DataOutput
        def taskOutputsResult = store.load("${hash}/outputs") as TaskOutputs
        then:
        taskRunResult == taskDescription
        dataOutputResult1 == dataOutput1
        dataOutputResult2 == dataOutput2
        taskOutputsResult.taskRun == "cid://1234567890"
        taskOutputsResult.workflowRun == "cid://hash"
        taskOutputsResult.outputs.size() == 3
        taskOutputsResult.outputs.get(0).type == FileOutParam.simpleName
        taskOutputsResult.outputs.get(0).name == "file1"
        taskOutputsResult.outputs.get(0).value == "cid://1234567890/outputs/fileOut1.txt"
        taskOutputsResult.outputs.get(1).type == FileOutParam.simpleName
        taskOutputsResult.outputs.get(1).name == "file2"
        taskOutputsResult.outputs.get(1).value == ["cid://1234567890/outputs/fileOut2.txt"]
        taskOutputsResult.outputs.get(2).type == ValueOutParam.simpleName
        taskOutputsResult.outputs.get(2).name == "id"
        taskOutputsResult.outputs.get(2).value == "value"

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task data output' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultCidStore();
        def session = Mock(Session) {
            getConfig()>>config
        }
        store.open(DataConfig.create(session))
        def observer = Spy(new CidObserver(session, store))
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
        def output = new DataOutput(outFile.toString(), new Checksum(fileHash, "nextflow", "standard"),
            "cid://15cd5b07", "cid://15cd5b07", attrs.size(), CidUtils.toDate(attrs.creationTime()), CidUtils.toDate(attrs.lastModifiedTime()) )
        and:
        observer.readAttributes(outFile) >> attrs

        when:
        observer.storeTaskOutput(task, outFile)
        then:
        folder.resolve(".meta/${hash}/outputs/foo/bar/file.bam/.data.json").text == new CidEncoder().encode(output)

        cleanup:
        folder?.deleteDir()
    }

    def 'should relativise task output dirs' (){
        when:
        def config = [workflow:[data:[enabled: true, store:[location:'cid']]]]
        def store = new DefaultCidStore();
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
        store.open(DataConfig.create(session))
        def observer = new CidObserver(session, store)
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
            def config = [workflow:[data:[enabled: true, store:[location:'cid']]]]
            def store = new DefaultCidStore();
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
        store.open(DataConfig.create(session))
        def observer = new CidObserver(session, store)
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
            def config = [workflow:[data:[enabled: true, store:[location:'cid']]]]
            def store = new DefaultCidStore();
            def session = Mock(Session) {
                getOutputDir()>>OUTPUT_DIR
                getConfig()>>config
            }
            store.open(DataConfig.create(session))
            def observer = new CidObserver(session, store)
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
            def config = [workflow:[data:[enabled: true, store:[location:'cid']]]]
            def store = new DefaultCidStore();
            def session = Mock(Session) {
                getOutputDir()>>OUTPUT_DIR
                getConfig()>>config
            }
            def observer = new CidObserver(session, store)
            observer.getWorkflowRelative(PATH)
        then:
            def e = thrown(IllegalArgumentException)
            e.message == "Cannot access relative path for workflow output '$PATH'"
        where:
        OUTPUT_DIR                      | PATH                                  | EXPECTED
        Path.of('/path/to/outDir')      | Path.of('/another/path/')             | "relative"
        Path.of('/path/to/outDir')      | Path.of('../relative')                | "relative"
    }

    def 'should save workflow output'() {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultCidStore();
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
        store.open(DataConfig.create(session))
        def observer = new CidObserver(session, store)
        def encoder = new CidEncoder()

        when: 'Starting workflow'
            observer.onFlowCreate(session)
            observer.onFlowBegin()
        then: 'History file should contain execution hash'
            def cid = store.getHistoryLog().getRecord(uniqueId).runCid.substring(CID_PROT.size())
            cid == observer.executionHash

        when: ' publish output with source file'
            def outFile1 = outputDir.resolve('foo/file.bam')
            Files.createDirectories(outFile1.parent)
            outFile1.text = 'some data1'
            def sourceFile1 = workDir.resolve('12/3987/file.bam')
            Files.createDirectories(sourceFile1.parent)
            sourceFile1.text = 'some data1'
            observer.onFilePublish(outFile1, sourceFile1)
            observer.onWorkflowPublish("a", outFile1)

        then: 'check file 1 output metadata in cid store'
            def attrs1 = Files.readAttributes(outFile1, BasicFileAttributes)
            def fileHash1 = CacheHelper.hasher(outFile1).hash().toString()
            def output1 = new DataOutput(outFile1.toString(), new Checksum(fileHash1, "nextflow", "standard"),
                "cid://123987/outputs/file.bam", "$CID_PROT${observer.executionHash}",
                attrs1.size(), CidUtils.toDate(attrs1.creationTime()), CidUtils.toDate(attrs1.lastModifiedTime()) )
            folder.resolve(".meta/${observer.executionHash}/outputs/foo/file.bam/.data.json").text == encoder.encode(output1)

        when: 'publish without source path'
        def outFile2 = outputDir.resolve('foo/file2.bam')
            Files.createDirectories(outFile2.parent)
            outFile2.text = 'some data2'
            def attrs2 = Files.readAttributes(outFile2, BasicFileAttributes)
            def fileHash2 = CacheHelper.hasher(outFile2).hash().toString()
            observer.onFilePublish(outFile2)
            observer.onWorkflowPublish("b", outFile2)
        then: 'Check outFile2 metadata in cid store'
            def output2 = new DataOutput(outFile2.toString(), new Checksum(fileHash2, "nextflow", "standard"),
                "cid://${observer.executionHash}" , "cid://${observer.executionHash}",
                attrs2.size(), CidUtils.toDate(attrs2.creationTime()), CidUtils.toDate(attrs2.lastModifiedTime()) )
            folder.resolve(".meta/${observer.executionHash}/outputs/foo/file2.bam/.data.json").text == encoder.encode(output2)

        when: 'Workflow complete'
            observer.onFlowComplete()
        then: 'Check history file is updated and Workflow Result is written in the cid store'
            def finalCid = store.getHistoryLog().getRecord(uniqueId).runCid.substring(CID_PROT.size())
            def resultsRetrieved = store.load("${finalCid}/outputs") as WorkflowOutputs
            resultsRetrieved.outputs == [a: "cid://${observer.executionHash}/outputs/foo/file.bam", b: "cid://${observer.executionHash}/outputs/foo/file2.bam"]

        cleanup:
            folder?.deleteDir()
    }

}
