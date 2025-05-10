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

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx

import static nextflow.lineage.fs.LinPath.*
/**
 * File to store a history of the workflow executions and their corresponding LIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class DefaultLinHistoryLog implements LinHistoryLog {

    Path path

    DefaultLinHistoryLog(Path folder) {
        this.path = folder
        if( !path.exists() )
            Files.createDirectories(path)
    }

    @Override
    void write(String name, UUID id, String launchLid, Date date = null) {
        assert name
        assert id
        final timestamp = date ?: new Date()
        final recordFile = path.resolve(launchId.substring(LID_PROT.size()))
        try {
            recordFile.text = new LinHistoryRecord(timestamp, name, id, null, launchLid, null).toString()
            log.trace("Record for $launchLid written in lineage history log ${FilesEx.toUriString(this.path)}")
        }catch (Throwable e) {
            log.warn("Can't write record $launchLid file ${FilesEx.toUriString(recordFile)}", e.message)
        }
    }

    @Override
    void finalize(String launchLid, String runLid, String status) {
        assert id
        final recordFile = path.resolve(launchId.substring(LID_PROT.size()))
        try {
            final current = LinHistoryRecord.parse(recordFile.text)
            recordFile.text = new LinHistoryRecord(current.timestamp, current.runName, id, status, current.launchLid, runLid).toString()
        }
        catch (Throwable e) {
            log.warn("Can't read record $launchId file: ${FilesEx.toUriString(recordFile)}", e.message)
        }
    }

    @Override
    List<LinHistoryRecord> getRecords(){
        List<LinHistoryRecord> list = new LinkedList<LinHistoryRecord>()
        try {
            this.path.eachFile { Path file -> list.add(LinHistoryRecord.parse(file.text))}
        }
        catch (Throwable e) {
            log.warn "Exception reading records from lineage history folder: ${FilesEx.toUriString(this.path)}", e.message
        }
        return list.sort {it.timestamp }
    }

}
