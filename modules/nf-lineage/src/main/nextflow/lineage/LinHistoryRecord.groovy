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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

import java.text.DateFormat
import java.text.SimpleDateFormat

/**
 * Record of workflow executions and their corresponding Lineage IDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'runName,sessionId')
class LinHistoryRecord {

    static final public DateFormat TIMESTAMP_FMT = new SimpleDateFormat('yyyy-MM-dd HH:mm:ss z')

    final Date timestamp
    final String runName
    final UUID sessionId
    final String runLid

    LinHistoryRecord(Date timestamp, String name, UUID sessionId, String runLid) {
        this.timestamp = timestamp
        this.runName = name
        this.sessionId = sessionId
        this.runLid = runLid
    }

    protected LinHistoryRecord() {}

    List<String> toList() {
        return List.of(
            timestamp ? TIMESTAMP_FMT.format(timestamp) : '-',
            runName ?: '-',
            sessionId.toString(),
            runLid ?: '-',
        )
    }

    @Override
    String toString() {
        toList().join('\t')
    }

    static LinHistoryRecord parse(String line) {
        final cols = line.tokenize('\t')
        if (cols.size() == 4) {
            return new LinHistoryRecord(TIMESTAMP_FMT.parse(cols[0]), cols[1], UUID.fromString(cols[2]), cols[3])
        }
        throw new IllegalArgumentException("Not a valid history entry: `$line`")
    }
}
