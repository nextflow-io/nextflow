package nextflow.sraql


import com.google.cloud.bigquery.BigQuery;
import com.google.cloud.bigquery.BigQueryOptions;
import com.google.cloud.bigquery.FieldValueList;
import com.google.cloud.bigquery.Job;
import com.google.cloud.bigquery.JobId;
import com.google.cloud.bigquery.JobInfo;
import com.google.cloud.bigquery.QueryJobConfiguration;
import com.google.cloud.bigquery.TableResult
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sraql.config.SraqlDataSource
import spock.lang.Specification


class Scratch extends Specification {

    def 'should perform query'() {
        given:
        def bigquery = BigQueryOptions.getDefaultInstance().getService();
        def queryConfig = QueryJobConfiguration.newBuilder(
                "SELECT * "
                        + " FROM `nih-sra-datastore.sra.metadata` as s "
                        + " WHERE organism = 'Mycobacterium tuberculosis'"
                        + " AND consent='public' "
                        + " LIMIT 10"
        )
                .setUseLegacySql(false)
                .build();

        def jobId = JobId.of(UUID.randomUUID().toString());
        def queryJob = bigquery
                .create(JobInfo.newBuilder(queryConfig).setJobId(jobId).build())
                .waitFor();


        when:
        def result = queryJob.getQueryResults();

        then:
        for (row in result.iterateAll()) {
            println(row)
        }

    }

}
