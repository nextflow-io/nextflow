package nextflow.sraql

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sraql.config.SraqlDataSource
import spock.lang.Specification
/**
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class QueryHandlerTest extends Specification {

    def 'should normalise query'() {
        given:
        def ext = new QueryHandler()

        expect:

        ext.normalize("SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis';") == "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis';"
        ext.normalize("SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis'") == "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis';"
    }


    def 'should perform the query'() {

        given:
        def dataSource = new SraqlDataSource(['source':'google-bigquery'])

        when:
        def result = new DataflowQueue()
        def query = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE  bioproject='PRJNA494931';"
        new QueryHandler()
                .withTarget(result)
                .withDataSource(dataSource)
                .withStatement(query)
                .perform()

        then:
        result.val.acc == 'SRR7974377'
        result.val.acc == 'SRR7974375'
        result.val.acc == 'SRR7974376'
        result.val == Channel.STOP

    }
}
