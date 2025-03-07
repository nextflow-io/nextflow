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

import groovy.transform.EqualsAndHashCode
import groovy.util.logging.Slf4j
import nextflow.util.WithLockFile

import java.nio.file.Path
import java.text.DateFormat
import java.text.SimpleDateFormat

/**
 * File to store a history of the workflow executions and their corresponding CIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
class CidHistoryFile extends WithLockFile {
    private static final DateFormat TIMESTAMP_FMT = new SimpleDateFormat('yyyy-MM-dd HH:mm:ss')

    CidHistoryFile(Path file) {
        super(file.toString())
    }

    void write(String name, UUID key, String runCid, Date date = null) {
        assert key

        withFileLock {
            def timestamp = date ?: new Date()
            log.debug("Writting record for $key in CID history file $this")
            this << new CidRecord(timestamp: timestamp, runName: name, sessionId: key, runCid: runCid).toString() << '\n'
        }
    }

    void update(UUID sessionId, String runCid) {
        assert sessionId

        try {
            withFileLock { update0(sessionId, runCid) }
        }
        catch (Throwable e) {
            log.warn "Can't update cid history file: $this", e
        }
    }

    String getRunCid(UUID id){
        assert id

        for (String line: this.readLines()){
            def current = line ? CidRecord.parse(line) : null
            if (current.sessionId == id) {
                return current.runCid
            }
        }
        log.warn("Can't find session $id in CID history file $this")
        return null
    }

    private void update0(UUID id, String runCid) {
        assert id
        def newHistory = new StringBuilder()

        this.readLines().each { line ->
            try {
                def current = line ? CidRecord.parse(line) : null
                if (current.sessionId == id) {
                    log.debug("Updating record for $id in CID history file $this")
                    current.runCid = runCid
                    newHistory << current.toString() << '\n'
                } else {
                    newHistory << line << '\n'
                }
            }
            catch (IllegalArgumentException e) {
                log.warn("Can't read CID history file: $this", e)
            }
        }

        // rewrite the history content
        this.setText(newHistory.toString())
    }

    @EqualsAndHashCode(includes = 'runName,sessionId')
    static class CidRecord {
        Date timestamp
        String runName
        UUID sessionId
        String runCid

        CidRecord(UUID sessionId, String name = null) {
            this.runName = name
            this.sessionId = sessionId
        }

        protected CidRecord() {}

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

        static CidRecord parse(String line) {
            def cols = line.tokenize('\t')
            if (cols.size() == 2)
                return new CidRecord(UUID.fromString(cols[0]))

            if (cols.size() == 4) {

                return new CidRecord(
                    timestamp: TIMESTAMP_FMT.parse(cols[0]),
                    runName: cols[1],
                    sessionId: UUID.fromString(cols[2]),
                    runCid: cols[3]
                )
            }

            throw new IllegalArgumentException("Not a valid history entry: `$line`")
        }
    }

}