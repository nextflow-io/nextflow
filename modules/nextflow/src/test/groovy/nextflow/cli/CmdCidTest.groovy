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

import groovy.json.JsonOutput

import java.nio.file.Files

import nextflow.SysEnv
import nextflow.dag.MermaidHtmlRenderer
import nextflow.data.cid.CidHistoryRecord
import nextflow.data.cid.CidStoreFactory
import nextflow.plugin.Plugins

import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 * CLI cid Tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CmdCidTest extends Specification {

    def setup() {
        // clear the environment to avoid the local env pollute the test env
        SysEnv.push([:])
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
            def folder = Files.createTempDirectory('test').toAbsolutePath()
            def configFile = folder.resolve('nextflow.config')
            configFile.text = "workflow.data.enabled = true\nworkflow.data.store.location = '$folder'".toString()
            def historyFile = folder.resolve(".meta/.history")
            Files.createDirectories(historyFile.parent)
            def uniqueId = UUID.randomUUID()
            def date = new Date();
            def launcher = Mock(Launcher){
                getOptions() >> new CliOptions(config: [configFile.toString()])
            }
            def recordEntry = "${CidHistoryRecord.TIMESTAMP_FMT.format(date)}\trun_name\t${uniqueId}\tcid://123456\tcid://456789".toString()
            historyFile.text = recordEntry
        when:
            def cidCmd = new CmdCid(launcher: launcher, args: ["log"])
            cidCmd.run()
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
        configFile.text = "workflow.data.enabled = true\nworkflow.data.store.location = '$folder'".toString()
        def historyFile = folder.resolve(".meta/.history")
        Files.createDirectories(historyFile.parent)
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        when:
        def cidCmd = new CmdCid(launcher: launcher, args: ["log"])
        cidCmd.run()
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

        cleanup:
        folder?.deleteDir()
    }

    def 'should show cid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "workflow.data.enabled = true\nworkflow.data.store.location = '$folder'".toString()
        def cidFile = folder.resolve(".meta/12345/.data.json")
        Files.createDirectories(cidFile.parent)
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }

        def recordEntry = JsonOutput.prettyPrint('{"type":"WorkflowOutput",' +
                '"path":"/path/to/file",' +
                '"checksum":"45372qe",' +
                '"source":"cid://123987/file.bam",' +
                '"size": 1234,' +
                '"createdAt": 123456789,' +
                '"modifiedAt": 123456789,' +
                '"annotations":null}')
        cidFile.text = recordEntry
        when:
            def cidCmd = new CmdCid(launcher: launcher, args: ["show", "cid://12345"])
            cidCmd.run()
            def stdout = capture
                .toString()
                .readLines()// remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .findResults { line -> !line.contains('INFO') ? line : null }
                .findResults { line -> !line.contains('plugin') ? line : null }

        then:
            stdout.size() == recordEntry.readLines().size()
            stdout.join('\n') == recordEntry

        cleanup:
            folder?.deleteDir()
    }

    def 'should warn if no cid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        configFile.text = "workflow.data.enabled = true\nworkflow.data.store.location = '$folder'".toString()
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }

        when:
            def cidCmd = new CmdCid(launcher: launcher, args: ["show", "cid://12345"])
            cidCmd.run()
            def stdout = capture
                .toString()
                .readLines()// remove the log part
                .findResults { line -> !line.contains('DEBUG') ? line : null }
                .findResults { line -> !line.contains('INFO') ? line : null }
                .findResults { line -> !line.contains('plugin') ? line : null }

        then:
            stdout.size() == 1
            stdout[0] == "No entry found for cid://12345."

        cleanup:
            folder?.deleteDir()
    }

    def 'should get lineage cid content' (){
        given:
        def folder = Files.createTempDirectory('test').toAbsolutePath()
        def configFile = folder.resolve('nextflow.config')
        def outputHtml = folder.resolve('lineage.html')
        configFile.text = "workflow.data.enabled = true\nworkflow.data.store.location = '$folder'".toString()
        def launcher = Mock(Launcher){
            getOptions() >> new CliOptions(config: [configFile.toString()])
        }
        def cidFile = folder.resolve(".meta/12345/file.bam/.data.json")
        def cidFile2 = folder.resolve(".meta/123987/file.bam/.data.json")
        def cidFile3 = folder.resolve(".meta/123987/.data.json")
        def cidFile4 = folder.resolve(".meta/45678/output.txt/.data.json")
        def cidFile5 = folder.resolve(".meta/45678/.data.json")
        Files.createDirectories(cidFile.parent)
        Files.createDirectories(cidFile2.parent)
        Files.createDirectories(cidFile3.parent)
        Files.createDirectories(cidFile4.parent)
        Files.createDirectories(cidFile5.parent)

        def recordEntry = JsonOutput.prettyPrint('{"type":"WorkflowOutput",' +
            '"path":"/path/to/file","checksum":"45372qe","source":"cid://123987/file.bam",' +
            '"size": 1234,"createdAt": 123456789, "modifiedAt": 123456789,"annotations":null}')
        cidFile.text = recordEntry
        recordEntry = JsonOutput.prettyPrint('{"type":"TaskOutput",' +
            '"path":"/path/to/file","checksum":"45372qe","source":"cid://123987",' +
            '"size": 1234,"createdAt": 123456789,"modifiedAt": 123456789,"annotations":null}')
        cidFile2.text = recordEntry
        recordEntry = JsonOutput.prettyPrint('{"type":"TaskRun",' +
            '"sessionId":"u345-2346-1stw2", "name":"foo","code":"abcde2345",' +
            '"inputs": [{"type": "ValueInParam","name": "sample_id","value": "ggal_gut"},' +
            '{"type": "FileInParam","name": "reads","value": ["cid://45678/output.txt"]}],' +
            '"container": null,"conda": null,"spack": null,"architecture": null,' +
            '"globalVars": {},"binEntries": [],"annotations":null}')
        cidFile3.text = recordEntry
        recordEntry = JsonOutput.prettyPrint('{"type":"TaskOutput",' +
            '"path":"/path/to/file","checksum":"45372qe","source":"cid://45678",' +
            '"size": 1234,"createdAt": 123456789,"modifiedAt": 123456789,"annotations":null}')
        cidFile4.text = recordEntry
        recordEntry = JsonOutput.prettyPrint('{"type":"TaskRun",' +
            '"sessionId":"u345-2346-1stw2", "name":"bar","code":"abfs2556",' +
            '"inputs": null,"container": null,"conda": null,"spack": null,"architecture": null,' +
            '"globalVars": {},"binEntries": [],"annotations":null}')
        cidFile5.text = recordEntry
        final network = """flowchart BT
    cid://12345/file.bam@{shape: document, label: "cid://12345/file.bam"}
    cid://123987/file.bam@{shape: document, label: "cid://123987/file.bam"}
    cid://123987@{shape: process, label: "foo"}
    ggal_gut@{shape: document, label: "ggal_gut"}
    cid://45678/output.txt@{shape: document, label: "cid://45678/output.txt"}
    cid://45678@{shape: process, label: "bar"}

    cid://123987/file.bam -->cid://12345/file.bam
    cid://123987 -->cid://123987/file.bam
    ggal_gut -->cid://123987
    cid://45678/output.txt -->cid://123987
    cid://45678 -->cid://45678/output.txt
"""
        final template = MermaidHtmlRenderer.readTemplate()
        def expectedOutput = template.replace('REPLACE_WITH_NETWORK_DATA', network)

        when:
        def cidCmd = new CmdCid(launcher: launcher, args: ["lineage", "cid://12345/file.bam", outputHtml.toString()])
        cidCmd.run()
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

        cleanup:
        folder?.deleteDir()

    }

}
