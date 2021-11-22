package nextflow.sraql


import com.google.cloud.bigquery.BigQueryOptions
import com.google.cloud.bigquery.Field
import com.google.cloud.bigquery.FieldValue;
import com.google.cloud.bigquery.JobId;
import com.google.cloud.bigquery.JobInfo;
import com.google.cloud.bigquery.QueryJobConfiguration;
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sraql.config.SraqlDataSource
import spock.lang.Specification


class Scratch extends Specification {

    def bigquery = BigQueryOptions.getDefaultInstance().getService()
    def queryConfig = QueryJobConfiguration.newBuilder(
            "SELECT * "
                    + " FROM `nih-sra-datastore.sra.metadata` as s "
                    + " WHERE s.organism = 'Mycobacterium tuberculosis'"
                    + " AND s.consent='public' "
                    + " AND s.sra_study='ERP124850' "
                    + " LIMIT 3"
    )
            .setUseLegacySql(false)
            .build()

    def jobId = JobId.of(UUID.randomUUID().toString());
    def queryJob = bigquery
            .create(JobInfo.newBuilder(queryConfig).setJobId(jobId).build())
            .waitFor()


    def 'should perform query'() {
        when:
        def result = queryJob.getQueryResults()

        then:

        for (row in result.iterateAll()) {
            def item = new ArrayList()

            for (field in result.schema.fields) {
                def fieldName = field.name
                def fieldContent = row.get(fieldName)

                if( fieldContent.attribute == FieldValue.Attribute.PRIMITIVE ) {
//                        println("${fieldName}: ${fieldContent.value}")
                    item <<  fieldContent.value
                }
            }
            println(item)
        }
    }
}
