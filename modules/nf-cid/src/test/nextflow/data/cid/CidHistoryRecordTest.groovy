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

/**
 * CID History Record tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CidHistoryRecordTest extends Specification {
    def "CidRecord parse should throw for invalid record"() {
        when:
        CidHistoryRecord.parse("invalid-record")

        then:
        thrown(IllegalArgumentException)
    }

    def "CidRecord parse should handle 4-column record"() {
        given:
        def timestamp = new Date()
        def formattedTimestamp = CidHistoryRecord.TIMESTAMP_FMT.format(timestamp)
        def line = "${formattedTimestamp}\trun-1\t${UUID.randomUUID()}\tcid://123"

        when:
        def record = CidHistoryRecord.parse(line)

        then:
        record.timestamp != null
        record.runName == "run-1"
        record.runCid == "cid://123"
    }

    def "CidRecord toString should produce tab-separated format"() {
        given:
        UUID sessionId = UUID.randomUUID()
        def record = new CidHistoryRecord(new Date(), "TestRun", sessionId, "cid://123")

        when:
        def line = record.toString()

        then:
        line.contains("\t")
        line.split("\t").size() == 4
    }
}
