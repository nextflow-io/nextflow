/*
 * Copyright 2020-2022, Seqera Labs
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
import java.sql.PreparedStatement
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
class QueryHandler implements QueryOp<QueryHandler> {

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
    private boolean emitColumns = false
    private Integer batchSize
    private long batchDelayMillis = 100
    private int queryCount

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

    QueryOp withOpts(Map opts) {
        if( opts.emitColumns )
            this.emitColumns = opts.emitColumns as boolean
        if( opts.batchSize )
            this.batchSize = opts.batchSize as Integer
        if( opts.batchDelay )
            this.batchDelayMillis = opts.batchDelay as long
        return this
    }

    int batchSize() {
        return batchSize
    }

    int queryCount() {
        return queryCount
    }

    @Override
    QueryHandler perform(boolean async=false) {
        final conn = connect(dataSource ?: SqlDataSource.DEFAULT)
        if( async )
            queryAsync(conn)
        else
            queryExec(conn)
        return this
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
        def future = CompletableFuture.runAsync ({ queryExec(conn) })
        future.exceptionally(this.&handlerException)
    }

    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }

    protected void queryExec(Connection conn) {
        if( batchSize ) {
            query1(conn)
        }
        else {
            query0(conn)
        }
    }

    protected void query0(Connection conn) {
        try {
            try (Statement stm = conn.createStatement()) {
                try( def rs = stm.executeQuery(normalize(statement)) ) {
                    if( emitColumns )
                        emitColumns(rs)
                    emitRowsAndClose(rs)
                }
            }
        }
        finally {
            conn.close()
        }
    }

    protected void query1(Connection conn) {
        try {
            // create the query adding the `offset` and `limit` params
            final query = makePaginationStm(statement)
            // create the prepared statement
            try (PreparedStatement stm = conn.prepareStatement(query)) {
                int count = 0
                int len = 0
                do {
                    final offset = (count++) * batchSize
                    final limit = batchSize

                    stm.setInt(1, limit)
                    stm.setInt(2, offset)
                    queryCount++
                    try ( def rs = stm.executeQuery() ) {
                        if( emitColumns && count==1 )
                            emitColumns(rs)
                        len = emitRows(rs)
                        sleep(batchDelayMillis)
                    }
                }
                while( len==batchSize )
            }
            finally {
                // close the channel
                target.bind(Channel.STOP)
            }
        }
        finally {
            conn.close()
        }
    }

    protected String makePaginationStm(String sql) {
        if( sql.toUpperCase().contains('LIMIT') )
                throw new IllegalArgumentException("Sql query should not include the LIMIT statement when pageSize is specified: $sql")
        if( sql.toUpperCase().contains('OFFSET') )
            throw new IllegalArgumentException("Sql query should not include the OFFSET statement when pageSize is specified: $sql")

        return sql.stripEnd(' ;') + " LIMIT ? OFFSET ?;"
    }

    protected emitColumns(ResultSet rs) {
        final meta = rs.getMetaData()
        final cols = meta.getColumnCount()

        def item = new ArrayList(cols)
        for( int i=0; i<cols; i++) {
            item[i] = meta.getColumnName(i+1)
        }
        // emit the value
        target.bind(item)
    }

    protected int emitRows(ResultSet rs) {
        final meta = rs.getMetaData()
        final cols = meta.getColumnCount()

        int count=0
        while( rs.next() ) {
            count++
            def item = new ArrayList(cols)
            for( int i=0; i<cols; i++) {
                item[i] = rs.getObject(i+1)
            }
            // emit the value
            target.bind(item)
        }
        return count
    }

    protected int emitRowsAndClose(ResultSet rs) {
        try {
            emitRows(rs)
        }
        finally {
            // close the channel
            target.bind(Channel.STOP)
        }
    }
}
