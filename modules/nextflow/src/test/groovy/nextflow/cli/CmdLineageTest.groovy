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

package nextflow.cli

import nextflow.SysEnv
import nextflow.dag.MermaidHtmlRenderer
import nextflow.lineage.DefaultLinHistoryLog
import nextflow.lineage.LinHistoryRecord
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.FileOutput
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.TaskRun
import nextflow.lineage.serde.LinEncoder
import nextflow.plugin.Plugins
import java.nio.file.Files
import java.time.OffsetDateTime
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 * CLI lineage Tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CmdLineageTest extends Specification {

    def setup() {
        // clear the environment to avoid the local env pollute the test env
        SysEnv.push([:])
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
            def folder = Files.createTempDirectory('test').toAbsolutePath()
            def configFile = folder.resolve('nextflow.config')
            configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
            def historyFile = folder.resolve(".history")
            def lidLog = new DefaultLinHistoryLog(historyFile)
            def uniqueId = UUID.randomUUID()
            def date = new Date();
            def launcher = Mock(Launcher){
                getOptions() >> new CliOptions(config: [configFile.toString()])
            }
            lidLog.write("run_name", uniqueId, "lid://123456", date)
            def recordEntry = "${LinHistoryRecord.TIMESTAMP_FMT.format(date)}\trun_name\t${uniqueId}\tlid://123456".toString()
        when:
            def lidCmd = new CmdLineage(launcher: launcher, args: ["list"])
            lidCmd.run()
            def stdout = capture
                .toString()
                .readLines()// remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .findResults { line -> !line.contains('INFO') ? line : null }
                .findResults { line -> !line.contains('plugin') ? line : null }

        then:
            stdout.size() == 2
            stdout[1] == recordEntry

        cleanup:
            folder?.deleteDir()
    }

    def 'should print no history' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
        def historyFile = folder.resolve(".history")
        Files.createDirectories(historyFile.parent)
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        when:
        def lidCmd = new CmdLineage(launcher: launcher, args: ["list"])
        lidCmd.run()
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

        cleanup:
        folder?.deleteDir()
    }

    def 'should show lid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
        def lidFile = folder.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        def time = OffsetDateTime.now()
        def encoder = new LinEncoder().withPrettyPrint(true)
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
                "lid://123987/file.bam","lid://12345/","lid://123987/", 1234, time, time, null)
        def jsonSer = encoder.encode(entry)
        def expectedOutput = jsonSer
        lidFile.text = jsonSer
        when:
            def lidCmd = new CmdLineage(launcher: launcher, args: ["view", "lid://12345"])
            lidCmd.run()
            def stdout = capture
                .toString()
                .readLines()// remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .findResults { line -> !line.contains('INFO') ? line : null }
                .findResults { line -> !line.contains('plugin') ? line : null }

        then:
            stdout.size() == expectedOutput.readLines().size()
            stdout.join('\n') == expectedOutput

        cleanup:
            folder?.deleteDir()
    }

    def 'should warn if no lid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }

        when:
            def lidCmd = new CmdLineage(launcher: launcher, args: ["view", "lid://12345"])
            lidCmd.run()
            def stdout = capture
                .toString()
                .readLines()// remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .findResults { line -> !line.contains('INFO') ? line : null }
                .findResults { line -> !line.contains('plugin') ? line : null }

        then:
            stdout.size() == 1
            stdout[0] == "Error loading lid://12345 - Lineage record 12345 not found"

        cleanup:
            folder?.deleteDir()
    }

    def 'should get lineage lid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        def outputHtml = folder.resolve('lineage.html')
        configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        def lidFile = folder.resolve("12345/file.bam/.data.json")
        def lidFile2 = folder.resolve("123987/file.bam/.data.json")
        def lidFile3 = folder.resolve("123987/.data.json")
        def lidFile4 = folder.resolve("45678/output.txt/.data.json")
        def lidFile5 = folder.resolve("45678/.data.json")
        Files.createDirectories(lidFile.parent)
        Files.createDirectories(lidFile2.parent)
        Files.createDirectories(lidFile3.parent)
        Files.createDirectories(lidFile4.parent)
        Files.createDirectories(lidFile5.parent)
        def encoder = new LinEncoder()
        def time = OffsetDateTime.now()
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
                "lid://123987/file.bam", "lid://45678",null, 1234, time, time, null)
        lidFile.text = encoder.encode(entry)
        entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
                "lid://123987", "lid://45678", "lid://123987", 1234, time, time, null)
        lidFile2.text = encoder.encode(entry)
        entry = new TaskRun("u345-2346-1stw2", "foo",
                new Checksum("abcde2345","nextflow","standard"),
                'this is a script',
                [new Parameter( "val", "sample_id","ggal_gut"),
                new Parameter("path","reads",["lid://45678/output.txt"])],
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
        final network = """\
            flowchart TB
                lid://12345/file.bam["lid://12345/file.bam"]
                lid://123987/file.bam["lid://123987/file.bam"]
                lid://123987(["foo [lid://123987]"])
                ggal_gut["ggal_gut"]
                lid://45678/output.txt["lid://45678/output.txt"]
                lid://45678(["bar [lid://45678]"])
                lid://123987/file.bam --> lid://12345/file.bam
                lid://123987 --> lid://123987/file.bam
                ggal_gut --> lid://123987
                lid://45678/output.txt --> lid://123987
                lid://45678 --> lid://45678/output.txt""".stripIndent()
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        def lidCmd = new CmdLineage(launcher: launcher, args: ["render", "lid://12345/file.bam", outputHtml.toString()])
        lidCmd.run()
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == 1
        stdout[0] == "Rendered lineage graph for lid://12345/file.bam to ${outputHtml}"
        outputHtml.exists()
        outputHtml.text == expectedOutput

        cleanup:
        folder?.deleteDir()

    }

    def 'should show query results'(){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "lineage.enabled = true\nlineage.store.location = '$folder'".toString()
        def lidFile = folder.resolve("12345/.data.json")
        Files.createDirectories(lidFile.parent)
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        def encoder = new LinEncoder().withPrettyPrint(true)
        def time = OffsetDateTime.now()
        def entry = new FileOutput("path/to/file",new Checksum("45372qe","nextflow","standard"),
                "lid://123987/file.bam", "lid://12345", "lid://123987/", 1234, time, time, ['foo', 'bar'])
        def jsonSer = encoder.encode(entry)
        def expectedOutput = '[\n  "lid://12345"\n]'
        lidFile.text = jsonSer
        when:
        def lidCmd = new CmdLineage(launcher: launcher, args: ["find", "type=FileOutput", "label=foo"])
        lidCmd.run()
        def stdout = capture
            .toString()
            .readLines()// remove the log part
            .findResults { line -> !line.contains('DEBUG') ? line : null }
            .findResults { line -> !line.contains('INFO') ? line : null }
            .findResults { line -> !line.contains('plugin') ? line : null }

        then:
        stdout.size() == expectedOutput.readLines().size()
        stdout.join('\n') == expectedOutput

        cleanup:
        folder?.deleteDir()
    }

}
