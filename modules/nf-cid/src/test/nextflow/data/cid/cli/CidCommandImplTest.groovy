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

package nextflow.data.cid.cli

import nextflow.SysEnv
import nextflow.config.ConfigMap
import nextflow.dag.MermaidHtmlRenderer
import nextflow.data.cid.CidHistoryRecord
import nextflow.data.cid.CidStoreFactory
import nextflow.data.cid.DefaultCidHistoryLog
import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.TaskRun
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.cid.serde.CidEncoder
import nextflow.plugin.Plugins
import org.junit.Rule
import spock.lang.Specification
import spock.lang.TempDir
import test.OutputCapture
import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant

class CidCommandImplTest extends Specification{

    @TempDir
    Path tmpDir

    Path storeLocation
    ConfigMap configMap

    def setup() {
        // clear the environment to avoid the local env pollute the test env
        SysEnv.push([:])
        storeLocation = tmpDir.resolve("store")
        configMap = new ConfigMap([workflow:[ data: [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]]])
    }

    def cleanup() {
        Plugins.stop()
        CidStoreFactory.reset()
        SysEnv.pop()
    }

    def setupSpec() {
        CidStoreFactory.reset()
    }
    /*
     * Read more http://mrhaki.blogspot.com.es/2015/02/spocklight-capture-and-assert-system.html
     */
    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should print executions cids' (){
        given:
        def historyFile = storeLocation.resolve(".meta/.history")
        def cidLog = new DefaultCidHistoryLog(historyFile)
        def uniqueId = UUID.randomUUID()
        def date = new Date();
        def recordEntry = "${CidHistoryRecord.TIMESTAMP_FMT.format(date)}\trun_name\t${uniqueId}\tcid://123456".toString()
        cidLog.write("run_name", uniqueId, "cid://123456", date)
        when:
        new CidCommandImpl().log(configMap)
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
        new CidCommandImpl().log(configMap)
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('WARN') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "No workflow runs CIDs found."

    }

    def 'should show cid content' (){
        given:
        def cidFile = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def time = Instant.ofEpochMilli(123456789)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam","cid://123987/", null, 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = jsonSer
        cidFile.text = jsonSer
        when:
        new CidCommandImpl().show(configMap, ["cid://12345"])
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

    def 'should warn if no cid content' (){
        given:

        when:
        new CidCommandImpl().show(configMap, ["cid://12345"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "No entries found for cid://12345"
    }

    def 'should get lineage cid content' (){
        given:

        def outputHtml = tmpDir.resolve('lineage.html')

        def cidFile = storeLocation.resolve(".meta/12345/file.bam/.data.json")
        def cidFile2 = storeLocation.resolve(".meta/123987/file.bam/.data.json")
        def cidFile3 = storeLocation.resolve(".meta/123987/.data.json")
        def cidFile4 = storeLocation.resolve(".meta/45678/output.txt/.data.json")
        def cidFile5 = storeLocation.resolve(".meta/45678/.data.json")
        Files.createDirectories(cidFile.parent)
        Files.createDirectories(cidFile2.parent)
        Files.createDirectories(cidFile3.parent)
        Files.createDirectories(cidFile4.parent)
        Files.createDirectories(cidFile5.parent)
        def encoder = new CidEncoder()
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://45678", null, 1234, time, time, null)
        cidFile.text = encoder.encode(entry)
        entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987", "cid://45678", "cid://123987", 1234, time, time, null)
        cidFile2.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "foo",
            new Checksum("abcde2345","nextflow","standard"),
            new Checksum("abfsc2375","nextflow","standard"),
            [new Parameter( "ValueInParam", "sample_id","ggal_gut"),
             new Parameter("FileInParam","reads",["cid://45678/output.txt"]),
             new Parameter("FileInParam","input",[new DataPath("path/to/file",new Checksum("45372qe","nextflow","standard"))])
            ],
            null, null, null, null, [:],[], null)
        cidFile3.text = encoder.encode(entry)
        entry  = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://45678", "cid://45678", null, 1234, time, time, null)
        cidFile4.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "bar",
            new Checksum("abfs2556","nextflow","standard"),
            new Checksum("abfsc2375","nextflow","standard"),
            null,null, null, null, null, [:],[], null)
        cidFile5.text = encoder.encode(entry)
        final network = """flowchart BT
    cid://12345/file.bam@{shape: document, label: "cid://12345/file.bam"}
    cid://123987/file.bam@{shape: document, label: "cid://123987/file.bam"}
    cid://123987@{shape: process, label: "foo"}
    ggal_gut@{shape: document, label: "ggal_gut"}
    path/to/file@{shape: document, label: "path/to/file"}
    cid://45678/output.txt@{shape: document, label: "cid://45678/output.txt"}
    cid://45678@{shape: process, label: "bar"}

    cid://123987/file.bam -->cid://12345/file.bam
    cid://123987 -->cid://123987/file.bam
    ggal_gut -->cid://123987
    cid://45678/output.txt -->cid://123987
    path/to/file -->cid://123987
    cid://45678 -->cid://45678/output.txt
"""
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new CidCommandImpl().lineage(configMap, ["cid://12345/file.bam", outputHtml.toString()])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Linage graph for cid://12345/file.bam rendered in ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput
    }

    def 'should get lineage from workflow cid content' (){
        given:

        def outputHtml = tmpDir.resolve('lineage.html')

        def cidFile = storeLocation.resolve(".meta/12345/file.bam/.data.json")
        def cidFile3 = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        Files.createDirectories(cidFile3.parent)
        def encoder = new CidEncoder()
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://12345", "cid://12345", null, 1234, time, time, null)
        cidFile.text = encoder.encode(entry)
        def wf = new Workflow([new DataPath("/path/to/main.nf)")], "hello-nf", "aasdklk")
        entry = new WorkflowRun(wf,"sessionId","run_name",
            [new Parameter( "String", "sample_id","ggal_gut"),
             new Parameter("Integer","reads",2)])
        cidFile3.text = encoder.encode(entry)
        final network = """flowchart BT
    cid://12345/file.bam@{shape: document, label: "cid://12345/file.bam"}
    cid://12345@{shape: processes, label: "run_name"}
    ggal_gut@{shape: document, label: "ggal_gut"}
    2.0@{shape: document, label: "2.0"}

    cid://12345 -->cid://12345/file.bam
    ggal_gut -->cid://12345
    2.0 -->cid://12345
"""
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        new CidCommandImpl().lineage(configMap, ["cid://12345/file.bam", outputHtml.toString()])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Linage graph for cid://12345/file.bam rendered in ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput
    }

    def 'should show query results'(){
        given:
        def cidFile = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://123987/", null, 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = jsonSer
        cidFile.text = jsonSer
        when:
        new CidCommandImpl().show(configMap, ["cid:///?type=DataOutput"])
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
        def cidFile = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def cidFile2 = storeLocation.resolve(".meta/67890/.data.json")
        Files.createDirectories(cidFile2.parent)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://123987/", null, 1234, time, time, null)
        def entry2 = new DataOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "cid://123987/file2.bam", "cid://123987/", null, 1235, time, time, null)
        def expectedOutput1 = '[\n  "path/to/file",\n  "path/to/file2"\n]'
        def expectedOutput2 = '[\n  "path/to/file2",\n  "path/to/file"\n]'
        cidFile.text = encoder.encode(entry)
        cidFile2.text = encoder.encode(entry2)
        when:
        new CidCommandImpl().show(configMap, ["cid:///?type=DataOutput#path"])
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
        def cidFile = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def cidFile2 = storeLocation.resolve(".meta/67890/.data.json")
        Files.createDirectories(cidFile2.parent)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://123987/", null, 1234, time, time, null)
        def entry2 = new DataOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "cid://123987/file2.bam", "cid://123987/", null, 1235, time, time, null)
        cidFile.text = encoder.encode(entry)
        cidFile2.text = encoder.encode(entry2)
        def expectedOutput = '''diff --git 12345 67890
--- 12345
+++ 67890
@@ -1,15 +1,15 @@
 {
   "type": "DataOutput",
-  "path": "path/to/file",
+  "path": "path/to/file2",
   "checksum": {
-    "value": "45372qe",
+    "value": "42472qet",
     "algorithm": "nextflow",
     "mode": "standard"
   },
-  "source": "cid://123987/file.bam",
+  "source": "cid://123987/file2.bam",
   "workflowRun": "cid://123987/",
   "taskRun": null,
-  "size": 1234,
+  "size": 1235,
   "createdAt": "1970-01-02T10:17:36.789Z",
   "modifiedAt": "1970-01-02T10:17:36.789Z",
   "annotations": null
'''

        when:
        new CidCommandImpl().diff(configMap, ["cid://12345", "cid://67890"])
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
        def cidFile = storeLocation.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://123987/", null, 1234, time, time, null)
        cidFile.text = encoder.encode(entry)

        when:
        new CidCommandImpl().diff(configMap, ["cid://89012", "cid://12345"])
        new CidCommandImpl().diff(configMap, ["cid://12345", "cid://67890"])
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 2
        stdout[0] == "No entry found for cid://89012."
        stdout[1] == "No entry found for cid://67890."
    }

    def 'should print error store is not found in diff'(){
        when:
        def config = new ConfigMap()
        new CidCommandImpl().log(config)
        new CidCommandImpl().show(config,["cid:///?type=DataOutput"])
        new CidCommandImpl().lineage(config,["cid://12345", "output.html"])
        new CidCommandImpl().diff(config, ["cid://89012", "cid://12345"])

        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }
        def expectedOutput = "Error CID store not loaded. Check Nextflow configuration."
        then:
        stdout.size() == 4
        stdout[0] == expectedOutput
        stdout[1] == expectedOutput
        stdout[2] == expectedOutput
        stdout[3] == expectedOutput
    }

    def 'should find metadata descriptions'(){
        given:
        def cidFile = storeLocation.resolve(".meta/123987/file.bam/.data.json")
        Files.createDirectories(cidFile.parent)
        def cidFile2 = storeLocation.resolve(".meta/123987/file2.bam/.data.json")
        Files.createDirectories(cidFile2.parent)
        def encoder = new CidEncoder().withPrettyPrint(true)
        def time = Instant.ofEpochMilli(123456789)
        def entry = new DataOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
            "cid://123987/file.bam", "cid://123987/", null, 1234, time, time, null)
        def entry2 = new DataOutput("path/to/file2",new Checksum("42472qet","nextflow","standard"),
            "cid://123987/file2.bam", "cid://123987/", null, 1235, time, time, null)
        def expectedOutput1 = '[\n  "cid://123987/file.bam",\n  "cid://123987/file2.bam"\n]'
        def expectedOutput2 = '[\n  "cid://123987/file2.bam",\n  "cid://123987/file.bam"\n]'
        cidFile.text = encoder.encode(entry)
        cidFile2.text = encoder.encode(entry2)
        when:
        new CidCommandImpl().find(configMap, ["type=DataOutput"])
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
