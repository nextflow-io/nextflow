/*
 * Copyright 2020-2021, Seqera Labs
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

package nextflow.sql

import java.sql.Connection
import java.sql.ResultSet
import java.sql.Statement
import java.util.concurrent.CompletableFuture

import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.sql.config.SqlDataSource

/**
 * Implement the logic for query a DB in async manner
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class QueryHandler implements QueryOp {

    private static Map<String,Class<?>> type_mapping = [:]

    static {
        type_mapping.CHAR = String
        type_mapping.VARCHAR = String
        type_mapping.LONGVARCHAR = String
        type_mapping.NUMERIC = BigDecimal
        type_mapping.DECIMAL = BigDecimal
        type_mapping.BIT = Boolean
        type_mapping.TINYINT = Byte
        type_mapping.SMALLINT = Short
        type_mapping.INTEGER = Integer
        type_mapping.BIGINT	= Long
        type_mapping.REAL= Float
        type_mapping.FLOAT= Double
        type_mapping.DOUBLE	= Double
        type_mapping.BINARY	= byte[]
        type_mapping.VARBINARY = byte[]
        type_mapping.LONGVARBINARY= byte[]
        type_mapping.DATE = java.sql.Date
        type_mapping.TIME = java.sql.Time
        type_mapping.TIMESTAMP= java.sql.Timestamp
    }

    private DataflowWriteChannel target
    private String statement
    private SqlDataSource dataSource

    @Override
    QueryOp withStatement(String stm) {
        this.statement = stm
        return this
    }

    @Override
    QueryOp withTarget(DataflowWriteChannel channel) {
        this.target = channel
        return this
    }

    @Override
    QueryOp withDataSource(SqlDataSource datasource) {
        this.dataSource = datasource
        return this
    }

    @Override
    void perform(boolean async=false) {
        final conn = connect(dataSource ?: SqlDataSource.DEFAULT)
        if( async )
            queryAsync(conn)
        else
            query0(conn)
    }

    protected Connection connect(SqlDataSource ds) {
        log.debug "Creating SQL connection: ${ds}"
        Sql.newInstance(ds.toMap()).getConnection()
    }

    protected String normalize(String q) {
        if( !q )
            throw new IllegalArgumentException("Missing query argument")
        def result = q.trim()
        if( !result.endsWith(';') )
            result += ';'
        return result
    }

    protected queryAsync(Connection conn) {
        def future = CompletableFuture.runAsync ({ query0(conn) })
        future.exceptionally(this.&handlerException)
    }

    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }

    protected void query0(Connection conn) {
        try {
            try (Statement stm = conn.createStatement()) {
                try( def rs = stm.executeQuery(normalize(statement)) ) {
                    emitRows0(rs)
                }
            }
        }
        finally {
            conn.close()
        }
    }

    protected emitRows0(ResultSet rs) {
        try {
            final meta = rs.getMetaData()
            final cols = meta.getColumnCount()

            while( rs.next() ) {
                def item = new ArrayList(cols)
                for( int i=0; i<cols; i++) {
                    item[i] = rs.getObject(i+1)
                }
                // emit the value
                target.bind(item)
            }
        }
        finally {
            // close the channel
            target.bind(Channel.STOP)
        }
    }
}
