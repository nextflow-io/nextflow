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

package nextflow.data.cid.h2

import java.sql.Timestamp

import com.zaxxer.hikari.HikariDataSource
import groovy.sql.GroovyRowResult
import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.data.cid.CidHistoryLog
import nextflow.data.cid.CidHistoryRecord

/**
 * Implement a {@link CidHistoryLog} based on H2 database
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class H2CidHistoryLog implements CidHistoryLog {

    private HikariDataSource dataSource

    H2CidHistoryLog(HikariDataSource dataSource) {
        this.dataSource = dataSource
    }

    @Override
    void write(String name, UUID sessionId, String runCid, String resultsCid) {
        try(final sql=new Sql(dataSource)) {
            def query = """
                INSERT INTO cid_history_record (timestamp, run_name, session_id, run_cid, results_cid) 
                VALUES (?, ?, ?, ?, ?)
            """
            def timestamp = new Timestamp(System.currentTimeMillis()) // Current timestamp
            sql.executeInsert(query, List.<Object>of(timestamp, name, sessionId.toString(), runCid, resultsCid))
        }
    }

    @Override
    void updateRunCid(UUID sessionId, String runCid) {
        try(final sql=new Sql(dataSource)) {
            def query = """
                UPDATE cid_history_record 
                SET run_cid = ? 
                WHERE session_id = ?
            """

            final count = sql.executeUpdate(query, List.<Object>of(runCid, sessionId.toString()))
            if (count > 0) {
                log.debug "Successfully updated run_cid for session_id: $sessionId"
            }
            else {
                log.warn "No record found with session_id: $sessionId"
            }
        }
    }

    @Override
    void updateResultsCid(UUID sessionId, String resultsCid) {
        try(final sql=new Sql(dataSource)) {
            def query = """
                UPDATE cid_history_record 
                SET results_cid = ? 
                WHERE session_id = ?
            """

            final count = sql.executeUpdate(query, List.<Object>of(resultsCid, sessionId.toString()))
            if (count > 0) {
                log.debug "Successfully updated run_cid for session_id: $sessionId"
            }
            else {
                log.warn "No record found with session_id: $sessionId"
            }
        }
    }

    @Override
    List<CidHistoryRecord> getRecords() {
        try(final sql=new Sql(dataSource)) {
            final result = new ArrayList<CidHistoryRecord>(100)
            final query = "SELECT * FROM cid_history_record "
            final rows = sql.rows(query)
            for( GroovyRowResult row : rows ) {
                result.add(
                    new CidHistoryRecord(
                        row.timestamp as Date,
                        row.run_name as String,
                        UUID.fromString(row.session_id as String),
                        row.run_cid as String,
                        row.results_cid as String
                    )
                )
            }
            return result
        }
    }

    @Override
    CidHistoryRecord getRecord(UUID sessionId) {
        try(final sql=new Sql(dataSource)) {
            final query = "SELECT * FROM cid_history_record WHERE session_id = ?"
            final row = sql.firstRow(query, sessionId.toString()) // Convert UUID to String for query
            if( !row )
                return null
            return new CidHistoryRecord(
                row.timestamp as Date,
                row.run_name as String,
                UUID.fromString(row.session_id as String),
                row.run_cid as String,
                row.results_cid as String
            )
        }
    }

}
