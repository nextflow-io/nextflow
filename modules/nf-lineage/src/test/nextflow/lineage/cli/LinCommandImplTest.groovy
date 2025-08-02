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

import nextflow.SysEnv
import nextflow.config.ConfigMap
import nextflow.dag.MermaidHtmlRenderer
import nextflow.lineage.LinHistoryRecord
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.DefaultLinHistoryLog
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.FileOutput
import nextflow.lineage.model.DataPath
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.TaskRun
import nextflow.lineage.model.Workflow
import nextflow.lineage.model.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.plugin.Plugins
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture
import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

class LinCommandImplTest extends Specification{

    @TempDir
    Path tmpDir

    Path storeLocation
    ConfigMap configMap

    def setup() {
        // clear the environment to avoid the local env pollute the test env
        SysEnv.push([:])
        storeLocation = tmpDir.resolve("store")
        configMap = new ConfigMap([lineage: [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]])
    }

    def cleanup() {
        Plugins.stop()
        LinStoreFactory.reset()
        SysEnv.pop()
    }

    def setupSpec() {
        LinStoreFactory.reset()
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
        new LinCommandImpl().log(configMap)
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 2
        stdout[1] == recordEntry
    }

    def 'should print no history' (){
        given:
        def historyFile = storeLocation.resolve(".meta/.history")
        Files.createDirectories(historyFile.parent)

        when:
        new LinCommandImpl().log(configMap)
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('WARN') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

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
        new LinCommandImpl().describe(configMap, ["lid://12345"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should warn if no lid content' (){
        given:

        when:
        new LinCommandImpl().describe(configMap, ["lid://12345"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Error loading lid://12345 - Lineage object 12345 not found"
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
            null, null, null, null, [:],[], null)
        lidFile3.text = encoder.encode(entry)
        entry  = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://45678", "lid://45678", null, 1234, time, time, null)
        lidFile4.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "bar",
            new Checksum("abfs2556","nextflow","standard"),
            'this is a script',
            null,null, null, null, null, [:],[], null)
        lidFile5.text = encoder.encode(entry)
        final network = """flowchart BT
    lid://12345/file.bam@{shape: document, label: "lid://12345/file.bam"}
    lid://123987/file.bam@{shape: document, label: "lid://123987/file.bam"}
    lid://123987@{shape: process, label: "foo [lid://123987]"}
    ggal_gut@{shape: document, label: "ggal_gut"}
    path/to/file@{shape: document, label: "path/to/file"}
    lid://45678/output.txt@{shape: document, label: "lid://45678/output.txt"}
    lid://45678@{shape: process, label: "bar [lid://45678]"}

    lid://123987/file.bam -->lid://12345/file.bam
    lid://123987 -->lid://123987/file.bam
    ggal_gut -->lid://123987
    lid://45678/output.txt -->lid://123987
    path/to/file -->lid://123987
    lid://45678 -->lid://45678/output.txt
"""
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new LinCommandImpl().render(configMap, ["lid://12345/file.bam", outputHtml.toString()])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Linage graph for lid://12345/file.bam rendered in ${outputHtml}"
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
        final network = """flowchart BT
    lid://12345/file.bam@{shape: document, label: "lid://12345/file.bam"}
    lid://12345@{shape: processes, label: "run_name [lid://12345]"}
    ggal_gut@{shape: document, label: "ggal_gut"}
    2.0@{shape: document, label: "2.0"}

    lid://12345 -->lid://12345/file.bam
    ggal_gut -->lid://12345
    2.0 -->lid://12345
"""
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new LinCommandImpl().render(configMap, ["lid://12345/file.bam", outputHtml.toString()])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Linage graph for lid://12345/file.bam rendered in ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput
    }

    def 'should show query results'(){
        given:
        def lidFile = storeLocation.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = jsonSer
        lidFile.text = jsonSer
        when:
        new LinCommandImpl().describe(configMap, ["lid:///?type=FileOutput"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput
    }

    def 'should show query with fragment'(){
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
        def expectedOutput1 = '[\n  "path/to/file",\n  "path/to/file2"\n]'
        def expectedOutput2 = '[\n  "path/to/file2",\n  "path/to/file"\n]'
        lidFile.text = encoder.encode(entry)
        lidFile2.text = encoder.encode(entry2)
        when:
        new LinCommandImpl().describe(configMap, ["lid:///?type=FileOutput#path"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.join('\n') == expectedOutput1 || stdout.join('\n') == expectedOutput2
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
@@ -1,15 +1,15 @@
 {
   "type": "FileOutput",
-  "path": "path/to/file",
+  "path": "path/to/file2",
   "checksum": {
-    "value": "45372qe",
+    "value": "42472qet",
     "algorithm": "nextflow",
     "mode": "standard"
   },
-  "source": "lid://123987/file.bam",
+  "source": "lid://123987/file2.bam",
   "workflowRun": "lid://123987/",
   "taskRun": null,
-  "size": 1234,
+  "size": 1235,
   "createdAt": "1970-01-02T10:17:36.789Z",
   "modifiedAt": "1970-01-02T10:17:36.789Z",
   "annotations": null
'''

        when:
        new LinCommandImpl().diff(configMap, ["lid://12345", "lid://67890"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

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
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 2
        stdout[0] == "No entry found for lid://89012."
        stdout[1] == "No entry found for lid://67890."
    }

    def 'should print error store is not found in diff'(){
        when:
        def config = new ConfigMap()
        new LinCommandImpl().log(config)
        new LinCommandImpl().describe(config, ["lid:///?type=FileOutput"])
        new LinCommandImpl().render(config, ["lid://12345", "output.html"])
        new LinCommandImpl().diff(config, ["lid://89012", "lid://12345"])

        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }
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
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(123456789), ZoneOffset.UTC)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "lid://123987/file.bam", "lid://123987/", null, 1234, time, time, null)
        def entry2 = new FileOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "lid://123987/file2.bam", "lid://123987/", null, 1235, time, time, null)
        def expectedOutput1 = '[\n  "lid://123987/file.bam",\n  "lid://123987/file2.bam"\n]'
        def expectedOutput2 = '[\n  "lid://123987/file2.bam",\n  "lid://123987/file.bam"\n]'
        lidFile.text = encoder.encode(entry)
        lidFile2.text = encoder.encode(entry2)
        when:
        new LinCommandImpl().find(configMap, ["type=FileOutput"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.join('\n') == expectedOutput1 || stdout.join('\n') == expectedOutput2
    }



}
