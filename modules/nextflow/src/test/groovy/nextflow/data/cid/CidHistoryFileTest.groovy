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
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path

/**
 * CID History file tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CidHistoryFileTest extends Specification {

    @TempDir
    Path tempDir

    Path historyFile
    CidHistoryFile cidHistoryFile

    def setup() {
        historyFile = tempDir.resolve("cid-history.txt")
        Files.createFile(historyFile)
        cidHistoryFile = new CidHistoryFile(historyFile)
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
        def parsedRecord = CidHistoryFile.CidRecord.parse(lines[0])
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
        def parsedRecord = CidHistoryFile.CidRecord.parse(lines[0])
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
        def parsedRecord = CidHistoryFile.CidRecord.parse(lines[0])
        parsedRecord.runCid == runCid
    }

    def "CidRecord parse should throw for invalid record"() {
        when:
        CidHistoryFile.CidRecord.parse("invalid-record")

        then:
        thrown(IllegalArgumentException)
    }

    def "CidRecord parse should handle 4-column record"() {
        given:
        def timestamp = new Date()
        def formattedTimestamp = CidHistoryFile.TIMESTAMP_FMT.format(timestamp)
        def line = "${formattedTimestamp}\trun-1\t${UUID.randomUUID()}\tcid://123"

        when:
        def record = CidHistoryFile.CidRecord.parse(line)

        then:
        record.timestamp != null
        record.runName == "run-1"
        record.runCid == "cid://123"
    }

    def "CidRecord toString should produce tab-separated format"() {
        given:
        UUID sessionId = UUID.randomUUID()
        def record = new CidHistoryFile.CidRecord(sessionId, "TestRun")
        record.timestamp = new Date()
        record.runCid = "cid://123"

        when:
        def line = record.toString()

        then:
        line.contains("\t")
        line.split("\t").size() == 4
    }
}

