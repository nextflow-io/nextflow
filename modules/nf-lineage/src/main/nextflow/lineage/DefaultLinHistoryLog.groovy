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

    void write(String name, UUID key, String runLid, Date date = null) {
        assert name
        assert key
        def timestamp = date ?: new Date()
        final recordFile = path.resolve(runLid.substring(LID_PROT.size()))
        try {
            recordFile.text = new LinHistoryRecord(timestamp, name, key, runLid).toString()
            log.trace("Record for $key written in lineage history log ${FilesEx.toUriString(this.path)}")
        }catch (Throwable e) {
            log.warn("Can't write record $key file ${FilesEx.toUriString(recordFile)}", e.message)
        }
    }

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
