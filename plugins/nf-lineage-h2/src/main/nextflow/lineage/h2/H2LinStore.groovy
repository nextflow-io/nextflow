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

import java.sql.Clob

import com.zaxxer.hikari.HikariDataSource
import groovy.json.JsonSlurper
import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.lineage.LinHistoryLog
import nextflow.lineage.LinStore
import nextflow.lineage.LinUtils
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.serde.LinSerializable
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.config.LineageStoreOpts
import nextflow.util.TestOnly
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class H2LinStore implements LinStore {

    private HikariDataSource dataSource
    private LinEncoder encoder

    @Override
    H2LinStore open(LineageConfig config) {
        assert config.store.location.startsWith('jdbc:h2:')
        log.info "Connecting CID H2 store: '${config.store.location}'"
        encoder = new LinEncoder()
        dataSource = createDataSource(config.store)
        // create the db tables
        createDbTables(dataSource)
        createAlias(dataSource)
        return this
    }

    static HikariDataSource createDataSource(LineageStoreOpts store) {
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
            CREATE TABLE IF NOT EXISTS lid_file (
                id BIGINT AUTO_INCREMENT PRIMARY KEY,  
                path VARCHAR UNIQUE NOT NULL, 
                metadata CLOB NOT NULL
            );
        
            CREATE TABLE IF NOT EXISTS lid_file_tag (
                file_id BIGINT NOT NULL, 
                tags TEXT NOT NULL,  
                PRIMARY KEY (file_id),
                FOREIGN KEY (file_id) REFERENCES lid_file(id) ON DELETE CASCADE
            );

            CREATE TABLE IF NOT EXISTS lid_history_record (
                id IDENTITY PRIMARY KEY,  -- Auto-increment primary key
                timestamp TIMESTAMP NOT NULL,
                run_name VARCHAR(255) NOT NULL,
                session_id UUID NOT NULL,
                run_lid VARCHAR(255) NOT NULL,
                UNIQUE (run_name, session_id) -- Enforce uniqueness constraint
            );
        ''')
        }
    }

    static void createAlias(HikariDataSource dataSource){
        try(final sql=new Sql(dataSource)) {
            sql.execute("""
            CREATE ALIAS IF NOT EXISTS JSON_MATCH FOR "nextflow.lineage.h2.H2LinStore.matchesJsonQuery"
            """)
        }
    }

    @Override
    void save(String key, LinSerializable object) {
        final value = encoder.encode(object)
        try(final sql=new Sql(dataSource)) {
            sql.execute("""
            INSERT INTO lid_file (path, metadata) VALUES (?, ?)
        """, [key, (Object)value])
        }
    }

    @Override
    LinSerializable load(String key) {
        try(final sql=new Sql(dataSource)) {
            final result = sql.firstRow("SELECT metadata FROM lid_file WHERE path = ?", List.<Object>of(key))
            return result ? encoder.decode(toValue(result.metadata).toString()) : null
        }
    }

    protected Object toValue(Object obj) {
        return obj instanceof Clob
            ? obj.characterStream.text
            : obj
    }

    @Override
    LinHistoryLog getHistoryLog() {
        return new H2LinHistoryLog(dataSource)
    }

    @Override
    Map<String, LinSerializable> search(String queryString) {
        final results= new HashMap<String, LinSerializable>()
        try(final sql=new Sql(dataSource)) {
            sql.eachRow("SELECT path, metadata FROM lid_file WHERE JSON_MATCH(metadata, ?)", List.<Object>of(queryString)) { row ->
                results.put(row['path'] as String, encoder.decode(toValue(row['metadata']) as String))
            }
        }
        return results
    }

    /**
     * JSON_MATCH implementation for h2
     * @param jsonString
     * @param queryString
     * @return
     */
    static boolean matchesJsonQuery(String jsonString, String queryString) {
        def json = new JsonSlurper().parseText(jsonString)
        def conditions = LinUtils.parseQuery(queryString)
        return LinUtils.checkParams(json, conditions)
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
