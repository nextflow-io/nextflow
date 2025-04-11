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

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx
/**
 * File to store a history of the workflow executions and their corresponding CIDs
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class DefaultCidHistoryLog implements CidHistoryLog {

    Path path

    DefaultCidHistoryLog(Path folder) {
        this.path = folder
        if( !path.exists() )
            Files.createDirectories(path)
    }

    void write(String name, UUID key, String runCid, Date date = null) {
        assert key
        def timestamp = date ?: new Date()
        final recordFile = path.resolve(key.toString())
        try {
            recordFile.text = new CidHistoryRecord(timestamp, name, key, runCid).toString()
            log.trace("Record for $key written in CID history log ${FilesEx.toUriString(this.path)}")
        }catch (Throwable e) {
            log.warn("Can't write record $key file ${FilesEx.toUriString(recordFile)}", e.message)
        }
    }

    void updateRunCid(UUID id, String runCid) {
        assert id
        final recordFile = path.resolve(id.toString())
        try {
            def current = CidHistoryRecord.parse(path.resolve(id.toString()).text)
            recordFile.text = new CidHistoryRecord(current.timestamp, current.runName, id, runCid).toString()
        }
        catch (Throwable e) {
            log.warn("Can't read session $id file: ${FilesEx.toUriString(recordFile)}", e.message)
        }
    }

    List<CidHistoryRecord> getRecords(){
        List<CidHistoryRecord> list = new LinkedList<CidHistoryRecord>()
        try {
            this.path.eachFile { Path file -> list.add(CidHistoryRecord.parse(file.text))}
        }
        catch (Throwable e) {
            log.warn "Exception reading records from CID history folder: ${FilesEx.toUriString(this.path)}", e.message
        }
        return list.sort {it.timestamp }
    }

    CidHistoryRecord getRecord(UUID id) {
        assert id
        final recordFile = path.resolve(id.toString())
        try {
            return CidHistoryRecord.parse(recordFile.text)
        } catch( Throwable e ) {
            log.warn("Can't find session $id in file: ${FilesEx.toUriString(recordFile)}", e.message)
            return null
        }
    }

}
