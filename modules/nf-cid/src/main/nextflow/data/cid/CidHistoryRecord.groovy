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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode

import java.text.DateFormat
import java.text.SimpleDateFormat

/**
 * Record of workflow executions and their corresponding CIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'runName,sessionId')
class CidHistoryRecord {
    public static final DateFormat TIMESTAMP_FMT = new SimpleDateFormat('yyyy-MM-dd HH:mm:ss')
    final Date timestamp
    final String runName
    final UUID sessionId
    final String runCid

    CidHistoryRecord(Date timestamp, String name, UUID sessionId, String runCid) {
        this.timestamp = timestamp
        this.runName = name
        this.sessionId = sessionId
        this.runCid = runCid
    }

    CidHistoryRecord(UUID sessionId, String name = null) {
        this.runName = name
        this.sessionId = sessionId
    }

    protected CidHistoryRecord() {}

    List<String> toList() {
        def line = new ArrayList<String>(4)
        line << (timestamp ? TIMESTAMP_FMT.format(timestamp) : '-')
        line << (runName ?: '-')
        line << (sessionId.toString())
        line << (runCid ?: '-')
    }

    @Override
    String toString() {
        toList().join('\t')
    }

    static CidHistoryRecord parse(String line) {
        def cols = line.tokenize('\t')
        if (cols.size() == 2)
            return new CidHistoryRecord(UUID.fromString(cols[0]))

        if (cols.size() == 4) {
            return new CidHistoryRecord(TIMESTAMP_FMT.parse(cols[0]), cols[1], UUID.fromString(cols[2]), cols[3])
        }

        throw new IllegalArgumentException("Not a valid history entry: `$line`")
    }
}
