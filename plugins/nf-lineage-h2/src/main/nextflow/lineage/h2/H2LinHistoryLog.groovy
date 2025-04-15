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
 *
 */

package nextflow.lineage.h2

import java.sql.Timestamp

import com.zaxxer.hikari.HikariDataSource
import groovy.sql.GroovyRowResult
import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.lineage.LinHistoryLog
import nextflow.lineage.LinHistoryRecord

/**
 * Implement a {@link LinHistoryLog} based on H2 database
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class H2LinHistoryLog implements LinHistoryLog {

    private HikariDataSource dataSource

    H2LinHistoryLog(HikariDataSource dataSource) {
        this.dataSource = dataSource
    }

    @Override
    void write(String name, UUID sessionId, String runCid) {
        try(final sql=new Sql(dataSource)) {
            def query = """
                INSERT INTO lid_history_record (timestamp, run_name, session_id, run_lid) 
                VALUES (?, ?, ?, ?)
            """
            def timestamp = new Timestamp(System.currentTimeMillis()) // Current timestamp
            sql.executeInsert(query, List.<Object>of(timestamp, name, sessionId.toString(), runCid))
        }
    }

    @Override
    void updateRunLid(UUID sessionId, String runLid) {
        try(final sql=new Sql(dataSource)) {
            def query = """
                UPDATE lid_history_record 
                SET run_lid = ? 
                WHERE session_id = ?
            """

            final count = sql.executeUpdate(query, List.<Object>of(runLid, sessionId.toString()))
            if (count > 0) {
                log.debug "Successfully updated run_lid for session_id: $sessionId"
            }
            else {
                log.warn "No record found with session_id: $sessionId"
            }
        }
    }

    @Override
    List<LinHistoryRecord> getRecords() {
        try(final sql=new Sql(dataSource)) {
            final result = new ArrayList<LinHistoryRecord>(100)
            final query = "SELECT * FROM lid_history_record "
            final rows = sql.rows(query)
            for( GroovyRowResult row : rows ) {
                result.add(
                    new LinHistoryRecord(
                        row.timestamp as Date,
                        row.run_name as String,
                        UUID.fromString(row.session_id as String),
                        row.run_lid as String,
                    )
                )
            }
            return result
        }
    }

    @Override
    LinHistoryRecord getRecord(UUID sessionId) {
        try(final sql=new Sql(dataSource)) {
            final query = "SELECT * FROM lid_history_record WHERE session_id = ?"
            final row = sql.firstRow(query, sessionId.toString()) // Convert UUID to String for query
            if( !row )
                return null
            return new LinHistoryRecord(
                row.timestamp as Date,
                row.run_name as String,
                UUID.fromString(row.session_id as String),
                row.run_lid as String,
            )
        }
    }

}
