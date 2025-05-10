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
package nextflow.lineage

import spock.lang.Specification

/**
 * Lineage History Record tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinHistoryRecordTest extends Specification {
    def "LinHistoryRecord parse should throw for invalid record"() {
        when:
        LinHistoryRecord.parse("invalid-record")

        then:
        thrown(IllegalArgumentException)
    }

    def "LinHistoryRecord parse should handle TSV record"() {
        given:
        def timestamp = new Date()
        def formattedTimestamp = LinHistoryRecord.TIMESTAMP_FMT.format(timestamp)
        def line = "${formattedTimestamp}\trun-1\t${UUID.randomUUID()}\tSUCCEEDED\tlid://123\tlid://456"

        when:
        def record = LinHistoryRecord.parse(line)

        then:
        record.timestamp != null
        record.runName == "run-1"
        record.status == "SUCCEEDED"
        record.launchLid == "lid://123"
        record.runLid == "lid://456"
    }

    def "LinHistoryRecord toString should produce tab-separated format"() {
        given:
        UUID sessionId = UUID.randomUUID()
        def record = new LinHistoryRecord(new Date(), "TestRun", sessionId, "SUCCEEDED", "lid://123", "lid://456")

        when:
        def line = record.toString()

        then:
        line.contains("\t")
        line.split("\t").size() == 6
    }
}
