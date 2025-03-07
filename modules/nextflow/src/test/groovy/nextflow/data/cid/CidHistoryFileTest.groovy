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
class CidHistoryFileTest extends Specification {

    Path tempDir
    Path historyFile
    CidHistoryFile cidHistoryFile

    def setup() {
        tempDir = Files.createTempDirectory("wdir")
        historyFile = tempDir.resolve("cid-history.txt")
        Files.createFile(historyFile)
        cidHistoryFile = new CidHistoryFile(historyFile)
    }

    def cleanup(){
        tempDir?.deleteDir()
    }

    def "write should append a new record to the file"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "TestRun"
        String runCid = "cid://123"

        when:
        cidHistoryFile.write(runName, sessionId, runCid)

        then:
        def lines = Files.readAllLines(historyFile)
        lines.size() == 1
        def parsedRecord = CidHistoryRecord.parse(lines[0])
        parsedRecord.sessionId == sessionId
        parsedRecord.runName == runName
        parsedRecord.runCid == runCid
    }

    def "getRunCid should return correct runCid for existing session"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"

        and:
        cidHistoryFile.write(runName, sessionId, runCid)

        expect:
        cidHistoryFile.getRunCid(sessionId) == runCid
    }

    def "getRunCid should return null if session does not exist"() {
        expect:
        cidHistoryFile.getRunCid(UUID.randomUUID()) == null
    }

    def "update should modify existing runCid for given session"() {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String initialCid = "cid-abc"
        String updatedCid = "cid-updated"

        and:
        cidHistoryFile.write(runName, sessionId, initialCid)

        when:
        cidHistoryFile.update(sessionId, updatedCid)

        then:
        def lines = Files.readAllLines(historyFile)
        lines.size() == 1
        def parsedRecord = CidHistoryRecord.parse(lines[0])
        parsedRecord.runCid == updatedCid
    }

    def "update should do nothing if session does not exist"() {
        given:
        UUID existingSessionId = UUID.randomUUID()
        UUID nonExistingSessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"

        and:
        cidHistoryFile.write(runName, existingSessionId, runCid)

        when:
        cidHistoryFile.update(nonExistingSessionId, "new-cid")

        then:
        def lines = Files.readAllLines(historyFile)
        lines.size() == 1
        def parsedRecord = CidHistoryRecord.parse(lines[0])
        parsedRecord.runCid == runCid
    }

    def 'should get records' () {
        given:
        UUID sessionId = UUID.randomUUID()
        String runName = "Run1"
        String runCid = "cid://123"
        and:
        cidHistoryFile.write(runName, sessionId, runCid)
        when:
        def records = cidHistoryFile.getRecords()
        then:
        records.size() == 1
        records[0].sessionId == sessionId
        records[0].runName == runName
        records[0].runCid == runCid
    }
}

