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

package nextflow.data.cid.h2


import java.nio.file.Path
import java.util.function.Consumer

import com.zaxxer.hikari.HikariDataSource
import groovy.sql.Sql
import groovy.transform.CompileStatic
import nextflow.data.cid.CidHistoryLog
import nextflow.data.cid.CidStore
import nextflow.data.config.DataConfig
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class H2CidStore implements CidStore {

    private HikariDataSource dataSource

    @Override
    void open(DataConfig config) {
        assert config.store.location.startsWith('jdbc:h2:')
        dataSource = new HikariDataSource()
        dataSource.jdbcUrl = config.store.location
        dataSource.username = 'sa'
        dataSource.password = ''
        dataSource.maximumPoolSize = 10
        // create DDL is missing
        try(final sql=new Sql(dataSource)) {
            sql.execute('''
            CREATE TABLE IF NOT EXISTS files (
                id BIGINT AUTO_INCREMENT PRIMARY KEY,  
                path VARCHAR UNIQUE NOT NULL, 
                metadata CLOB NOT NULL
            );
        
            CREATE TABLE IF NOT EXISTS file_tags (
                file_id BIGINT NOT NULL, 
                tags TEXT NOT NULL,  
                PRIMARY KEY (file_id),
                FOREIGN KEY (file_id) REFERENCES files(id) ON DELETE CASCADE
            );    
        ''')
        }
    }

    @Override
    void save(String key, Object value) {
        try(final sql=new Sql(dataSource)) {
            sql.execute("""
            INSERT INTO files (path, metadata) VALUES (?, ?)
        """, [key, value])
        }
    }

    @Override
    void list(String key, Consumer<String> consumer) {

    }

    @Override
    Object load(String key) {
        try(final sql=new Sql(dataSource)) {
            return sql.firstRow("SELECT * FROM files WHERE path = ?", List.of(key))
        }
    }

    @Override
    Path getPath() {
        return null
    }

    @Override
    CidHistoryLog getHistoryLog() {
        return null
    }

    @Override
    void close() {
        dataSource.close()
    }

}
