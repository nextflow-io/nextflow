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

import java.sql.Clob

import com.zaxxer.hikari.HikariDataSource
import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.data.cid.CidHistoryLog
import nextflow.data.cid.CidStore
import nextflow.data.cid.serde.CidEncoder
import nextflow.data.cid.serde.CidSerializable
import nextflow.data.config.DataConfig
import nextflow.data.config.DataStoreOpts
import nextflow.util.TestOnly
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class H2CidStore implements CidStore {

    private HikariDataSource dataSource
    private CidEncoder encoder

    @Override
    H2CidStore open(DataConfig config) {
        assert config.store.location.startsWith('jdbc:h2:')
        log.info "Connecting CID H2 store: '${config.store.location}'"
        encoder = new CidEncoder()
        dataSource = createDataSource(config.store)
        // create the db tables
        createDbTables(dataSource)
        return this
    }

    static HikariDataSource createDataSource(DataStoreOpts store) {
        final result = new HikariDataSource()
        result.jdbcUrl = store.location
        result.driverClassName = 'org.h2.Driver'
        result.username = 'sa'
        result.password = ''
        result.maximumPoolSize = 10
        return result
    }

    static void createDbTables(HikariDataSource dataSource) {
        // create DDL is missing
        try(final sql=new Sql(dataSource)) {
            sql.execute('''
            CREATE TABLE IF NOT EXISTS cid_file (
                id BIGINT AUTO_INCREMENT PRIMARY KEY,  
                path VARCHAR UNIQUE NOT NULL, 
                metadata CLOB NOT NULL
            );
        
            CREATE TABLE IF NOT EXISTS cid_file_tag (
                file_id BIGINT NOT NULL, 
                tags TEXT NOT NULL,  
                PRIMARY KEY (file_id),
                FOREIGN KEY (file_id) REFERENCES cid_file(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS cid_history_record (
                id IDENTITY PRIMARY KEY,  -- Auto-increment primary key
                timestamp TIMESTAMP NOT NULL,
                run_name VARCHAR(255) NOT NULL,
                session_id UUID NOT NULL,
                run_cid VARCHAR(255) NOT NULL,
                results_cid VARCHAR(255) NOT NULL,
                UNIQUE (run_name, session_id) -- Enforce uniqueness constraint
            );
        ''')
        }
    }

    @Override
    void save(String key, CidSerializable object) {
        final value = encoder.encode(object)
        try(final sql=new Sql(dataSource)) {
            sql.execute("""
            INSERT INTO cid_file (path, metadata) VALUES (?, ?)
        """, [key, (Object)value])
        }
    }

    @Override
    CidSerializable load(String key) {
        try(final sql=new Sql(dataSource)) {
            final result = sql.firstRow("SELECT metadata FROM cid_file WHERE path = ?", List.<Object>of(key))
            return result ? encoder.decode(toValue(result.metadata).toString()) : null
        }
    }

    protected Object toValue(Object obj) {
        return obj instanceof Clob
            ? obj.characterStream.text
            : obj
    }

    @Override
    CidHistoryLog getHistoryLog() {
        return new H2CidHistoryLog(dataSource)
    }

    @Override
    void close() {
        dataSource.close()
    }

    @TestOnly
    void truncateAllTables() {
        try(final sql=new Sql(dataSource)) {
            println "Truncating all tables..."
            sql.execute("SET REFERENTIAL_INTEGRITY FALSE") // Disable foreign key constraints

            def tables = sql.rows("SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'PUBLIC'")
            tables.each { table ->
                final stm = "TRUNCATE TABLE ${table.TABLE_NAME}" as String
                sql.execute(stm) // Truncate each table
            }

            sql.execute("SET REFERENTIAL_INTEGRITY TRUE") // Re-enable constraints
            println "All tables truncated successfully"
        }
    }
}
