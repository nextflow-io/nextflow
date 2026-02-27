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
 */

package nextflow.lineage.cli

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import nextflow.SysEnv
import nextflow.config.ConfigMap
import nextflow.dag.MermaidHtmlRenderer
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.lineage.DefaultLinHistoryLog
import nextflow.lineage.LinHistoryRecord
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.fs.LinFileSystemProvider
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.plugin.Plugins
import nextflow.util.CacheHelper
import org.junit.Rule
import spock.lang.Shared
import spock.lang.Specification
import test.OutputCapture

import static test.TestHelper.filterLogNoise

class LinCommandImplTest extends Specification{

    @Shared
    Path tmpDir

    @Shared
    Path storeLocation

    @Shared
    ConfigMap configMap

    def reset() {
        def provider = FileHelper.getProviderFor('lid') as LinFileSystemProvider
        provider?.reset()
        LinStoreFactory.reset()
    }

    def setup() {
        reset()
        // clear the environment to avoid the local env pollute the test env
        SysEnv.push([:])
        tmpDir = Files.createTempDirectory('tmp')
        storeLocation = tmpDir.resolve("store")
        configMap = new ConfigMap([lineage: [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]])
    }

    def cleanup() {
        Plugins.stop()
        LinStoreFactory.reset()
        SysEnv.pop()
        tmpDir?.deleteDir()
    }

    def setupSpec() {
        reset()
    }

    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should print executions lids' (){
        given:
        def historyFile = storeLocation.resolve(".history")
        def lidLog = new DefaultLinHistoryLog(historyFile)
        def uniqueId = UUID.randomUUID()
        def date = new Date();
        def recordEntry = "${LinHistoryRecord.TIMESTAMP_FMT.format(date)}\trun_name\t${uniqueId}\tlid://123456".toString()
        lidLog.write("run_name", uniqueId, "lid://123456", date)
        when:
        new LinCommandImpl().list(configMap)
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == 2
        stdout[1] == recordEntry
    }

    def 'should print no history' (){
        given:
        def historyFile = storeLocation.resolve(".meta/.history")
        Files.createDirectories(historyFile.parent)

        when:
        new LinCommandImpl().list(configMap)
        def stdout = filterLogNoise(capture, true)

        then:
        stdout.size() == 1
        stdout[0] == "No workflow runs found in lineage history log"
    }

    def 'should show lid content' (){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam","lid://123987/", null, 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = jsonSer
        lidFile.text = jsonSer
        when:
        new LinCommandImpl().view(configMap, ["lid://12345"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should show empty lists content' (){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam","lid://123987/", null, 1234, time, time, [])
        def jsonSer = encoder.encode(entry)
        def expectedOutput = '[]'
        lidFile.text = jsonSer
        when:
        new LinCommandImpl().view(configMap, ["lid://12345#labels"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should show empty lists when no outputs' () {
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def expectedOutput = '[]'
        def wf = new Workflow([new DataPath("/path/to/main.nf)")], "hello-nf", "aasdklk")
        def entry = new WorkflowRun(wf, "sessionId", "run_name",
            [new Parameter("String", "sample_id", "ggal_gut"),
             new Parameter("Integer", "reads", 2)])
        lidFile.text = encoder.encode(entry)

        when:
        new LinCommandImpl().view(configMap, ["lid://12345#output"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should warn if no lid content' (){
        given:

        when:
        new LinCommandImpl().view(configMap, ["lid://12345"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == 1
        stdout[0] == "Error loading lid://12345 - Lineage record 12345 not found"
    }

    def 'should get lineage lid content' (){
        given:

        def outputHtml = tmpDir.resolve('lineage.html')

        def lidFile = storeLocation.resolve("12345/file.bam/.data.json")
        def lidFile2 = storeLocation.resolve("123987/file.bam/.data.json")
        def lidFile3 = storeLocation.resolve("123987/.data.json")
        def lidFile4 = storeLocation.resolve("45678/output.txt/.data.json")
        def lidFile5 = storeLocation.resolve("45678/.data.json")
        Files.createDirectories(lidFile.parent)
        Files.createDirectories(lidFile2.parent)
        Files.createDirectories(lidFile3.parent)
        Files.createDirectories(lidFile4.parent)
        Files.createDirectories(lidFile5.parent)
        def encoder = new LinEncoder()
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://45678", null, 1234, time, time, null)
        lidFile.text = encoder.encode(entry)
        entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987", "lid://45678", "lid://123987", 1234, time, time, null)
        lidFile2.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "foo",
            new Checksum("abcde2345","nextflow","standard"),
            'this is a script',
            [new Parameter( "val", "sample_id","ggal_gut"),
             new Parameter("path","reads", ["lid://45678/output.txt"] ),
             new Parameter("path","input", [new DataPath("path/to/file",new Checksum("45372qe","nextflow","standard"))])
            ],
            null, null, null, null, [:],[])
        lidFile3.text = encoder.encode(entry)
        entry  = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://45678", "lid://45678", null, 1234, time, time, null)
        lidFile4.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "bar",
            new Checksum("abfs2556","nextflow","standard"),
            'this is a script',
            null,null, null, null, null, [:],[])
        lidFile5.text = encoder.encode(entry)
        final network = """\
            flowchart TB
                lid://12345/file.bam["lid://12345/file.bam"]
                lid://123987/file.bam["lid://123987/file.bam"]
                lid://123987(["foo [lid://123987]"])
                ggal_gut["ggal_gut"]
                path/to/file["path/to/file"]
                lid://45678/output.txt["lid://45678/output.txt"]
                lid://45678(["bar [lid://45678]"])
                lid://123987/file.bam --> lid://12345/file.bam
                lid://123987 --> lid://123987/file.bam
                ggal_gut --> lid://123987
                lid://45678/output.txt --> lid://123987
                path/to/file --> lid://123987
                lid://45678 --> lid://45678/output.txt""".stripIndent()
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new LinCommandImpl().render(configMap, ["lid://12345/file.bam", outputHtml.toString()])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == 1
        stdout[0] == "Rendered lineage graph for lid://12345/file.bam to ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput
    }

    def 'should get lineage from workflow lid content' (){
        given:

        def outputHtml = tmpDir.resolve('lineage.html')

        def lidFile = storeLocation.resolve("12345/file.bam/.data.json")
        def lidFile3 = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        Files.createDirectories(lidFile3.parent)
        def encoder = new LinEncoder()
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://12345", "lid://12345", null, 1234, time, time, null)
        lidFile.text = encoder.encode(entry)
        def wf = new Workflow([new DataPath("/path/to/main.nf)")], "hello-nf", "aasdklk")
        entry = new WorkflowRun(wf,"sessionId","run_name",
            [new Parameter( "String", "sample_id","ggal_gut"),
             new Parameter("Integer","reads",2)])
        lidFile3.text = encoder.encode(entry)
        final network = """\
            flowchart TB
                lid://12345/file.bam["lid://12345/file.bam"]
                lid://12345(["run_name [lid://12345]"])
                ggal_gut["ggal_gut"]
                2.0["2.0"]
                lid://12345 --> lid://12345/file.bam
                ggal_gut --> lid://12345
                2.0 --> lid://12345""".stripIndent()
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new LinCommandImpl().render(configMap, ["lid://12345/file.bam", outputHtml.toString()])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == 1
        stdout[0] == "Rendered lineage graph for lid://12345/file.bam to ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput
    }

    def 'should show an error if trying to do a query'(){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = "Error loading lid:///?type=FileOutput - Cannot get record from the root LID URI"
        lidFile.text = jsonSer
        when:
        new LinCommandImpl().view(configMap, ["lid:///?type=FileOutput"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should diff'(){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def lidFile2 = storeLocation.resolve("67890/.data.json")
        Files.createDirectories(lidFile2.parent)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, null)
        def entry2 = new FileOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "lid://123987/file2.bam", "lid://123987/", null, 1235, time, time, null)
        lidFile.text = encoder.encode(entry)
        lidFile2.text = encoder.encode(entry2)
        def expectedOutput = '''diff --git 12345 67890
--- 12345
+++ 67890
@@ -2,16 +2,16 @@
   "version": "lineage/v1beta1",
   "kind": "FileOutput",
   "spec": {
-    "path": "path/to/file",
+    "path": "path/to/file2",
     "checksum": {
-      "value": "45372qe",
+      "value": "42472qet",
       "algorithm": "nextflow",
       "mode": "standard"
     },
-    "source": "lid://123987/file.bam",
+    "source": "lid://123987/file2.bam",
     "workflowRun": "lid://123987/",
     "taskRun": null,
-    "size": 1234,
+    "size": 1235,
     "createdAt": "1970-01-02T10:17:36.789Z",
     "modifiedAt": "1970-01-02T10:17:36.789Z",
     "labels": null
'''

        when:
        new LinCommandImpl().diff(configMap, ["lid://12345", "lid://67890"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.join('\n') == expectedOutput
    }

    def 'should print error if no entry found diff'(){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, null)
        lidFile.text = encoder.encode(entry)

        when:
        new LinCommandImpl().diff(configMap, ["lid://89012", "lid://12345"])
        new LinCommandImpl().diff(configMap, ["lid://12345", "lid://67890"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.size() == 2
        stdout[0] == "No entry found for lid://89012."
        stdout[1] == "No entry found for lid://67890."
    }

    def 'should print error store is not found in diff'(){
        when:
        def config = new ConfigMap()
        new LinCommandImpl().list(config)
        new LinCommandImpl().view(config, ["lid:///12345"])
        new LinCommandImpl().render(config, ["lid://12345", "output.html"])
        new LinCommandImpl().diff(config, ["lid://89012", "lid://12345"])

        def stdout = filterLogNoise(capture)
        def expectedOutput = "Error lineage store not loaded - Check Nextflow configuration"
        then:
        stdout.size() == 4
        stdout[0] == expectedOutput
        stdout[1] == expectedOutput
        stdout[2] == expectedOutput
        stdout[3] == expectedOutput
    }

    def 'should find metadata descriptions'(){
        given:
        def lidFile = storeLocation.resolve("123987/file.bam/.data.json")
        Files.createDirectories(lidFile.parent)
        def lidFile2 = storeLocation.resolve("123987/file2.bam/.data.json")
        Files.createDirectories(lidFile2.parent)
        def lidFile3 = storeLocation.resolve(".meta/123987/file3.bam/.data.json")
        Files.createDirectories(lidFile3.parent)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, ["experiment=test"])
        def entry2 = new FileOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "lid://123987/file2.bam", "lid://123987/", null, 1235, time, time, ["experiment=test"])
        def entry3 = new FileOutput("path/to/file3",new Checksum("42472qet","nextflow","standard"),
            "lid://123987/file2.bam", "lid://123987/", null, 1235, time, time, null)
        def expectedOutput1 = 'lid://123987/file.bam\nlid://123987/file2.bam'
        def expectedOutput2 = 'lid://123987/file2.bam\nlid://123987/file.bam'
        lidFile.text = encoder.encode(entry)
        lidFile2.text = encoder.encode(entry2)
        lidFile3.text = encoder.encode(entry3)
        when:
        new LinCommandImpl().find(configMap, ["type=FileOutput", "label=experiment=test"])
        def stdout = filterLogNoise(capture)

        then:
        stdout.join('\n') == expectedOutput1 || stdout.join('\n') == expectedOutput2
    }

    def 'should print correct validate path' () {
        given:
        def outputFolder = tmpDir.resolve('output')
        Files.createDirectories(outputFolder)
        def outputFile = outputFolder.resolve('file1.txt')
        outputFile.text = "this   is file1  == "

        and:
        def encoder = new LinEncoder().withPrettyPrint(true)
        def hash = CacheHelper.hasher(outputFile).hash().toString()
        def correctData = new FileOutput(outputFile.toString(), new Checksum(hash,"nextflow", "standard"))
        def incorrectData = new FileOutput(outputFile.toString(), new Checksum("incorrectHash","nextflow", "standard"))
        def lid1 = storeLocation.resolve('12345/output/file1.txt/.data.json')
        Files.createDirectories(lid1.parent)
        lid1.text = encoder.encode(correctData)
        def lid2 = storeLocation.resolve('12345/output/file2.txt/.data.json')
        Files.createDirectories(lid2.parent)
        lid2.text = encoder.encode(incorrectData)
        when:
        new LinCommandImpl().check(configMap, ["lid://12345/output/file1.txt"])
        def stdout = filterLogNoise(capture)
        def expectedOutput1 = "Checksum validation succeed"
        then:
        stdout.size() == 1
        stdout[0] == expectedOutput1

        when:
        new LinCommandImpl().check(configMap, ["lid://12345/output/file2.txt"])
        then:
        def err = thrown(AbortOperationException)
        err.message == "Checksum of '${outputFile}' does not match with lineage metadata"
    }

}
