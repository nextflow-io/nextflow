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

import groovy.json.JsonOutput
import nextflow.data.config.DataConfig
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.script.ScriptBinding
import nextflow.script.WorkflowMetadata
import nextflow.util.CacheHelper
import nextflow.util.PathNormalizer

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import spock.lang.Specification

import static nextflow.data.cid.fs.CidPath.CID_PROT

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidObserverTest extends Specification {

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
        def expectedString = '{"type":"WorkflowRun","workflow":{"type": "Workflow",' +
            '"mainScriptFile":{"path":"file://' + scriptFile.toString() +
            '", "checksum": {"value": "78910", "algorithm": "nextflow", "mode": "standard"}},' +
            '"otherScriptFiles": [], "repository": "https://nextflow.io/nf-test/",' +
            '"commitId": "123456" },' +
            '"sessionId": "' + uniqueId + '",' +
            '"name": "test_run", "params": []}'
        when:
        observer.onFlowCreate(session)
        observer.onFlowBegin()
        then:
        folder.resolve(".meta/${observer.executionHash}/.data.json").text == JsonOutput.prettyPrint(expectedString)

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task run' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[enabled: true, store:[location:folder.toString()]]]]
        def store = new DefaultCidStore();
        def uniqueId = UUID.randomUUID()
        def session = Mock(Session) {
            getConfig()>>config
            getUniqueId()>>uniqueId
            getRunName()>>"test_run"
        }
        store.open(DataConfig.create(session))
        def observer = new CidObserver(session, store)
        and:
        def hash = HashCode.fromInt(123456789)
        and:
        def processor = Mock(TaskProcessor){
            getTaskGlobalVars(_) >> [:]
            getTaskBinEntries(_) >> []
        }
        def task = Mock(TaskRun) {
            getId() >> TaskId.of(100)
            getName() >> 'foo'
            getHash() >> hash
            getProcessor() >> processor
            getSource() >> 'echo task source'
        }
        def sourceHash =CacheHelper.hasher('echo task source').hash().toString()
        def normalizer = Mock(PathNormalizer.class) {
            normalizePath( _ as Path) >> {Path p -> p?.toString()}
            normalizePath( _ as String) >> {String p -> p}
        }
        def expectedString = '{"type":"TaskRun",' +
            '"sessionId":"'+uniqueId.toString() + '",' +
            '"name":"foo", "codeChecksum": {' +
            '"value": "' + sourceHash + '", "algorithm": "nextflow", "mode": "standard"},' +
            '"inputs": null,"container": null,"conda": null,' +
            '"spack": null,"architecture": null,' +
            '"globalVars": {},"binEntries": [],"annotations":null}'
        when:
        observer.storeTaskRun(task, normalizer)
        then:
        folder.resolve(".meta/${hash.toString()}/.data.json").text == JsonOutput.prettyPrint(expectedString)

        cleanup:
        folder?.deleteDir()
    }

    def 'should save task output' () {
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
        def expectedString = '{"type":"TaskOutput",' +
            '"path":"' + outFile.toString() + '",' +
            '"checksum": { "value":"'+ fileHash + '",' +
            '"algorithm": "nextflow", "mode": "standard"},' +
            '"source":"cid://15cd5b07",' +
            '"size":'+attrs.size() + ',' +
            '"createdAt":' + attrs.creationTime().toMillis() + ',' +
            '"modifiedAt":'+ attrs.lastModifiedTime().toMillis() + ',' +
            '"annotations":null}'

        and:
        observer.readAttributes(outFile) >> attrs

        when:
        observer.storeTaskOutput(task, outFile)
        then:
        folder.resolve(".meta/${hash}/foo/bar/file.bam/.data.json").text
            == JsonOutput.prettyPrint(expectedString)

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

    def 'should return exception when relativize task output dirs' (){
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
            def e = thrown(Exception)
            e.message == "Cannot asses the relative path for output $PATH of ${task.name}".toString()

        where:
        WORK_DIR                            | STORE_DIR                         | PATH
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('/another/path/relative')
        Path.of('/path/to/work/12/3456789') | Path.of('/path/to/storeDir')      | Path.of('../path/to/storeDir/relative')
    }

    def 'should relativise workflow output dirs' (){
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

    def 'should return exception when relativise workflow output dirs' (){
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
            def e = thrown(Exception)
            e.message == "Cannot asses the relative path for workflow output $PATH"
        where:
        OUTPUT_DIR                      | PATH                                  | EXPECTED
        Path.of('/path/to/outDir')      | Path.of('/another/path/')             | "relative"
        Path.of('/path/to/outDir')      | Path.of('../relative')                | "relative"


    }

    def 'should save workflow output' (){
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
        then: 'check file 1 output metadata in cid store'
            def attrs1 = Files.readAttributes(outFile1, BasicFileAttributes)
            def fileHash1 = CacheHelper.hasher(outFile1).hash().toString()
            def expectedString1 =  '{"type":"WorkflowOutput",' +
                '"path":"' + outFile1.toString() + '",' +
                '"checksum": {"value": "'+ fileHash1 + '",' +
                '"algorithm": "nextflow", "mode": "standard"},' +
                '"source":"cid://123987/file.bam",' +
                '"size":'+attrs1.size() + ',' +
                '"createdAt":' + attrs1.creationTime().toMillis() + ',' +
                '"modifiedAt":'+ attrs1.lastModifiedTime().toMillis() + ',' +
                '"annotations":null}'
            folder.resolve(".meta/${observer.executionHash}/foo/file.bam/.data.json").text == JsonOutput.prettyPrint(expectedString1)

        when: 'publish without source path'
        def outFile2 = outputDir.resolve('foo/file2.bam')
            Files.createDirectories(outFile2.parent)
            outFile2.text = 'some data2'
            def attrs2 = Files.readAttributes(outFile2, BasicFileAttributes)
            def fileHash2 = CacheHelper.hasher(outFile2).hash().toString()
            observer.onFilePublish(outFile2)
        then: 'Check outFile2 metadata in cid store'
            def expectedString2 =  '{"type":"WorkflowOutput",' +
                '"path":"' + outFile2.toString() + '",' +
                '"checksum": { "value": "'+ fileHash2 + '",' +
                '"algorithm": "nextflow", "mode": "standard"},' +
                '"source":"cid://' + observer.executionHash +'",' +
                '"size":'+attrs2.size() + ',' +
                '"createdAt":' + attrs2.creationTime().toMillis() + ',' +
                '"modifiedAt":'+ attrs2.lastModifiedTime().toMillis() + ',' +
                '"annotations":null}'
            folder.resolve(".meta/${observer.executionHash}/foo/file2.bam/.data.json").text == JsonOutput.prettyPrint(expectedString2)

        when: 'Workflow complete'
            observer.onFlowComplete()
        then: 'Check history file is updated and Workflow Result is written in the cid store'
            def expectedString3 =  '{"type":"WorkflowResults",' +
                '"run":"cid://' + observer.executionHash +'",' +
                '"outputs": {}, "publishedFiles": [ ' +
                '"cid://'+ observer.executionHash + '/foo/file.bam",' +
                '"cid://'+ observer.executionHash + '/foo/file2.bam" ]}'
            def finalCid = store.getHistoryLog().getRecord(uniqueId).resultsCid.substring(CID_PROT.size())
            finalCid != observer.executionHash
            folder.resolve(".meta/${finalCid}/.data.json").text == JsonOutput.prettyPrint(expectedString3)

        cleanup:
            folder?.deleteDir()
    }

}
