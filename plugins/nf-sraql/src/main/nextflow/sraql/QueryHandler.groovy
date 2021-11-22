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

package nextflow.sraql

import com.google.cloud.bigquery.BigQuery
import com.google.cloud.bigquery.BigQueryOptions
import com.google.cloud.bigquery.FieldValue
import com.google.cloud.bigquery.JobId
import com.google.cloud.bigquery.JobInfo
import com.google.cloud.bigquery.QueryJobConfiguration
import com.google.cloud.bigquery.TableResult

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
import nextflow.sraql.config.SraqlDataSource

/**
 * Implement the logic for query an SRAQL datasource in async manner
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@Slf4j
@CompileStatic
class QueryHandler implements QueryOp {

    private DataflowWriteChannel target
    private String statement
    private SraqlDataSource dataSource

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
    QueryOp withDataSource(SraqlDataSource datasource) {
        this.dataSource = datasource
        return this
    }

    @Override
    void perform(boolean async = false) {
        final conn = connect(dataSource)
        if( async )
            queryAsync(conn)
        else
            query0(conn)
    }

    protected BigQuery connect(SraqlDataSource ds) {
        log.debug "Creating SRAQL connection: ${ds.source}"
        BigQueryOptions.getDefaultInstance().getService();
    }

    protected String normalize(String q) {
        if( !q )
            throw new IllegalArgumentException("Missing query argument")
        def result = q.trim()
        if( !result.endsWith(';') )
            result += ';'
        return result
    }

    protected queryAsync(BigQuery conn) {
        def future = CompletableFuture.runAsync({ query0(conn) })
        future.exceptionally(this.&handlerException)
    }

    static private void handlerException(Throwable e) {
        final error = e.cause ?: e
        log.error(error.message, error)
        final session = Global.session as Session
        session?.abort(error)
    }

    protected void query0(BigQuery conn) {

        def jobId = JobId.of(UUID.randomUUID().toString());

        def queryJobConfig = QueryJobConfiguration.newBuilder(
                normalize(statement)
        )
                .setUseLegacySql(false)
                .build();

        def queryJob = conn
                .create(JobInfo.newBuilder(queryJobConfig).setJobId(jobId).build())
                .waitFor()

        emitRows0(queryJob.getQueryResults())
    }

    protected emitRows0(TableResult result) {
        try {

            for (row in result.iterateAll()) {
                def fieldValues = new ArrayList()

                for (field in result.schema.fields) {
                    def fieldName = field.name
                    def fieldContent = row.get(fieldName)

                    if( fieldContent.attribute == FieldValue.Attribute.PRIMITIVE ) {
                        fieldValues << [fieldName, fieldContent.value]
                    }
                }
                target.bind(fieldValues)
            }
        }
        finally {
            // close the channel
            target.bind(Channel.STOP)
        }
    }
}
