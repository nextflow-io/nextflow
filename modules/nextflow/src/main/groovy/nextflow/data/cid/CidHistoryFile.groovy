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

import groovy.util.logging.Slf4j

import java.nio.channels.FileLock
import java.nio.file.Path

/**
 * File to store a history of the workflow executions and their corresponding CIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
class CidHistoryFile implements CidHistoryLog {

    Path path

    CidHistoryFile(Path file) {
        log.debug("History file $file")
        this.path = file
    }

    void write(String name, UUID key, String runCid, Date date = null) {
        assert key

        withFileLock {
            def timestamp = date ?: new Date()
            log.debug("Writting record for $key in CID history file $this")
            path << new CidHistoryRecord(timestamp: timestamp, runName: name, sessionId: key, runCid: runCid).toString() << '\n'
        }
    }

    void update(UUID sessionId, String runCid) {
        assert sessionId

        try {
            withFileLock { update0(sessionId, runCid) }
        }
        catch (Throwable e) {
            log.warn "Can't update CID history file: $this", e.message
        }
    }

    List<CidHistoryRecord> getRecords(){
        List<CidHistoryRecord> list = new LinkedList<CidHistoryRecord>()
        try {
            withFileLock { this.path.eachLine {list.add(CidHistoryRecord.parse(it)) } }
        }
        catch (Throwable e) {
            log.warn "Can't read records from CID history file: $this", e.message
        }
        return list

    }


    String getRunCid(UUID id) {
        assert id

        for (String line : this.path.readLines()) {
            def current = line ? CidHistoryRecord.parse(line) : null
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

        this.path.readLines().each { line ->
            try {
                def current = line ? CidHistoryRecord.parse(line) : null
                if (current.sessionId == id) {
                    log.debug("Updating record for $id in CID history file $this")
                    current.runCid = runCid
                    newHistory << current.toString() << '\n'
                } else {
                    newHistory << line << '\n'
                }
            }
            catch (IllegalArgumentException e) {
                log.warn("Can't read CID history file: $this", e.message)
            }
        }

        // rewrite the history content
        this.path.setText(newHistory.toString())
    }

    /**
     * Apply the given action by using a file lock
     *
     * @param action The closure implementing the action to be executed with a file lock
     * @return The value returned by the action closure
     */
    protected withFileLock(Closure action) {

        def rnd = new Random()
        long ts = System.currentTimeMillis()
        final parent = this.path.parent ?: Path.of('.').toAbsolutePath()
        def file = parent.resolve("${this.path.name}.lock".toString()).toFile()
        def fos = new FileOutputStream(file)
        try {
            Throwable error
            FileLock lock = null

            try {
                while (true) {
                    lock = fos.getChannel().tryLock()
                    if (lock) break
                    if (System.currentTimeMillis() - ts < 1_000)
                        sleep rnd.nextInt(75)
                    else {
                        error = new IllegalStateException("Can't lock file: ${this.path.toAbsolutePath()} -- Nextflow needs to run in a file system that supports file locks")
                        break
                    }
                }
                if (lock) {
                    return action.call()
                }
            }
            catch (Exception e) {
                return action.call()
            }
            finally {
                if (lock?.isValid()) lock.release()
            }

            if (error) throw error
        }
        finally {
            fos.closeQuietly()
            file.delete()
        }
    }
}