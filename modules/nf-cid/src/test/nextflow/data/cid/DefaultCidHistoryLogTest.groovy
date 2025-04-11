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

import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Path

/**
 * CID History file tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultCidHistoryLogTest extends Specification {

    Path tempDir
    Path historyFile
    DefaultCidHistoryLog cidHistoryLog

    def setup() {
        tempDir = Files.createTempDirectory("wdir")
        historyFile = tempDir.resolve("cid-history")
        cidHistoryLog = new DefaultCidHistoryLog(historyFile)
    }

    def cleanup(){
        tempDir?.deleteDir()
    }

    def "write should add a new file to the history folder"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "TestRun"
        String runCid = "cid://123"

        when:
        cidHistoryLog.write(runName, sessionId, runCid)

        then:
        def files = historyFile.listFiles()
        files.size() == 1
        def parsedRecord = CidHistoryRecord.parse(files[0].text)
        parsedRecord.sessionId == sessionId
        parsedRecord.runName == runName
        parsedRecord.runCid == runCid
    }

    def "should return correct record for existing session"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"

        and:
        cidHistoryLog.write(runName, sessionId, runCid)

        when:
        def record = cidHistoryLog.getRecord(sessionId)
        then:
        record.sessionId == sessionId
        record.runName == runName
        record.runCid == runCid
    }

    def "should return null and warn if session does not exist"() {
        expect:
        cidHistoryLog.getRecord(UUID.randomUUID()) == null
    }

    def "update should modify existing Cid for given session"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCidUpdated = "run-cid-updated"
        String resultsCidUpdated = "results-cid-updated"

        and:
        cidHistoryLog.write(runName, sessionId, 'run-cid-initial')

        when:
        cidHistoryLog.updateRunCid(sessionId, runCidUpdated)

        then:
        def files = historyFile.listFiles()
        files.size() == 1
        def parsedRecord = CidHistoryRecord.parse(files[0].text)
        parsedRecord.runCid == runCidUpdated
    }

    def "update should do nothing if session does not exist"() {
        given:
        UUID existingSessionId = UUID.randomUUID()
        UUID nonExistingSessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"
        and:
        cidHistoryLog.write(runName, existingSessionId, runCid)

        when:
        cidHistoryLog.updateRunCid(nonExistingSessionId, "new-cid")
        then:
        def files = historyFile.listFiles()
        files.size() == 1
        def parsedRecord = CidHistoryRecord.parse(files[0].text)
        parsedRecord.runCid == runCid
    }

    def 'should get records' () {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"
        and:
        cidHistoryLog.write(runName, sessionId, runCid)

        when:
        def records = cidHistoryLog.getRecords()
        then:
        records.size() == 1
        records[0].sessionId == sessionId
        records[0].runName == runName
        records[0].runCid == runCid
    }
}

