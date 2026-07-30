/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.lineage

import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Path

/**
 * Lineage History file tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultLinHistoryLogTest extends Specification {

    Path tempDir
    Path historyFile
    DefaultLinHistoryLog linHistoryLog

    def setup() {
        tempDir = Files.createTempDirectory("wdir")
        historyFile = tempDir.resolve("lin-history")
        linHistoryLog = new DefaultLinHistoryLog(historyFile)
    }

    def cleanup(){
        tempDir?.deleteDir()
    }

    def "write should add a new file to the history folder"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "TestRun"
        String runLid = "lid://123"

        when:
        linHistoryLog.write(runName, sessionId, runLid)

        then:
        def files = historyFile.listFiles()
        files.size() == 1
        def parsedRecord = LinHistoryRecord.parse(files[0].text)
        parsedRecord.sessionId == sessionId
        parsedRecord.runName == runName
        parsedRecord.runLid == runLid
    }

    def 'should get records' () {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runLid = "lid://123"
        and:
        linHistoryLog.write(runName, sessionId, runLid)

        when:
        def records = linHistoryLog.getRecords()
        then:
        records.size() == 1
        records[0].sessionId == sessionId
        records[0].runName == runName
        records[0].runLid == runLid
    }
}

