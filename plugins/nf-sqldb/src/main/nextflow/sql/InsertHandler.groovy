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

import groovy.sql.Sql
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import groovy.util.logging.Slf4j
import nextflow.sql.config.SqlDataSource
import nextflow.util.TupleHelper
/**
 * Handle SQL insert operation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class InsertHandler implements Closeable {

    static final List<String> QM = ['?']
    static final int DEFAULT_BATCH_SIZE = 10

    private Map opts
    private SqlDataSource ds
    private Connection connection
    private String sqlStatement
    private String into
    private List<String> columns
    private int batchSize
    private int batchCount
    private PreparedStatement preparedStatement
    private String setupStatement

    InsertHandler(SqlDataSource ds, Map opts) {
        this.ds = ds
        this.opts = opts ?: Collections.emptyMap()
        this.into = this.opts.into
        this.columns = cols0(this.opts.columns)
        this.sqlStatement = this.opts.statement
        this.batchSize = this.opts.batch ? this.opts.batch as int : DEFAULT_BATCH_SIZE
        this.batchSize = this.opts.batchSize ? this.opts.batchSize as int : DEFAULT_BATCH_SIZE
        this.setupStatement = this.opts.setup
        if( this.opts.batch )
            log.warn "The option 'batch' for the 'sqlInsert' operator has been deprecated - Use 'batchSize' instead"
        if( batchSize<1 )
            throw new IllegalArgumentException("SQL batch option must be greater than zero: $batchSize")
    }

    private Connection getConnection() {
        if( connection == null ) {
            connection = Sql.newInstance(ds.toMap()).getConnection()
            checkCreate(connection)
            connection.setAutoCommit(false)
        }
        return connection
    }

    private void checkCreate(Connection connection) {
        if( !setupStatement )
            return
        // this allows performing a SQL statement ideally create the required table
        // note: the underlying DB should support the idiom `create table if not exist`
        try( def stm = connection.createStatement()) {
            stm.executeUpdate(setupStatement)
        }
    }

    protected List cols0(names) {
        if( !names )
            return null
        if( names instanceof List )
            return names
        if( names instanceof CharSequence )
            return names.tokenize(',').collect(it -> it.trim())
        throw new IllegalArgumentException("Invalid sqlInsert 'columns' option - Offending value: $names [${names.getClass().getName()}]")
    }

    protected String makeStatement(Object entry) {
        if( sqlStatement )
            return sqlStatement
        if( !into )
            throw new IllegalArgumentException("Missing sqlInsert target table - Use 'into' option to specify it")

        if( !columns )
            columns = fetchColNames(entry)
        if( !columns )
            throw new IllegalArgumentException("Missing sqlInsert 'columns' option")

        final cnt = columns.size()
        final result = "insert into $into ( ${columns.join(',')} ) values ( ${(QM * cnt).join(',')} )"
        return result
    }

    protected List<String> fetchColNames(Object entry) {
        if( entry instanceof List ) {
            return fetchColNames0(into)
        }
        else if( entry instanceof Map ) {
            return fetchColNames1(into, (Map)entry)
        }
        else
            throw new IllegalStateException("Invalid sqlInsert entry data type: $entry [${entry?.getClass()?.getName()}]")
    }

    protected String getSqlStatement(Object entry) {
        if( sqlStatement )
            return sqlStatement
        sqlStatement = makeStatement(entry)
        log.debug "[SQL] sqlInsert statement: \"${sqlStatement}\""
        return sqlStatement
    }

    void perform(Object entry) {
        // make sure is not null
        if( entry==null )
            throw new IllegalArgumentException("Invalid sqlInsert entry - It must be a non-null value")

        if( entry instanceof Map ) {
            final sql = getSqlStatement(entry)
            performAsMap(sql, entry)
        }
        else if( entry instanceof List ) {
            final sql = getSqlStatement(entry)
            performAsTuple(sql, entry)
        }
        else {
            entry = TupleHelper.listOf(entry)
            final sql = getSqlStatement(entry)
            performAsTuple(sql, entry)
        }

    }

    protected void performAsMap(String sql, Map record) {
        executeStm(sql) {stm ->
            // loop over the expected columns and fetch a corresponding record value to be inserted
            for(int i=0; i<columns.size(); i++ ) {
                final col = columns[i]
                final value = record.get(col)
                stm.setObject(i+1, value)
            }
            // report a debug line
            log.debug "[SQL] perform sql statemet=$sql; entry=$record"
        }
    }

    protected void performAsTuple(String sql, List tuple) {
        executeStm(sql) { stm ->
            // loop over the tuple values and set a corresponding sql statement value
            for(int i=0; i<tuple.size(); i++ ) {
                def value = tuple[i]
                stm.setObject(i+1, value)
            }
            // report a debug line
            log.debug "[SQL] perform sql statemet=$sql; entry=$tuple"
        }
    }

    protected void executeStm(String sql, @ClosureParams(value = SimpleType.class, options = "java.sql.PreparedStatement") Closure operation) {
        // create the statement if needed
        if( preparedStatement==null )
            preparedStatement = getConnection().prepareStatement(sql)

        // make sure to clear previous params
        preparedStatement.clearParameters()
        // invoke the operation closure with setup the statement params
        operation.call(preparedStatement)
        // add to current batch
        preparedStatement.addBatch()
        // if the batch size is reached, perform the batch statement
        if( ++batchCount == batchSize ) {
            try {
                preparedStatement.executeBatch()
                preparedStatement.clearBatch()
            }
            finally {
                // reset the current batch count
                batchCount = 0
                // make sure to commit the current batch
                connection.commit()
            }
        }
    }

    @Memoized
    protected List<String> fetchColNames0(String tableName) {
        final conn = getConnection()
        final metaData  = conn.getMetaData()
        final resultSet = metaData.getColumns(null, null, tableName, null);

        final result = [] as List<String>
        while (resultSet.next()){
            result.add( resultSet.getString("COLUMN_NAME") )
        }

        return result
    }

    protected List<String> fetchColNames1(String tableName, Map<String,Object> record) {
        final keys = record.keySet().collect(it -> it.toUpperCase())
        final allCols = fetchColNames0(tableName)
        final result = new ArrayList(allCols.size())
        // return only DB columns that belongs map keys
        for( String it : allCols ) {
            if( it in keys )
                result.add(it)
        }
        // check if the record contains fields not mapped to table columns
        final diff = keys - result
        if( diff )
            throw new IllegalArgumentException("Unknown columns for table '$into': ${diff.join(',')} - Offending record: $record")
        return result
    }

    void close() {
        try {
            if( preparedStatement && batchCount>0 ) {
                log.debug("[SQL] flushing and committing open batch")
                preparedStatement.executeBatch()
                preparedStatement.close()
                connection.commit()
            }
        }
        finally {
            if( connection ) {
                log.debug "[SQL] closing JDBC connection"
                connection.close()
            }
        }
    }
}
