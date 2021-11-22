package nextflow.sraql

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sraql.QueryHandler
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
        when:
        def result = new DataflowQueue()
        def query = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis'  AND sra_study='ERP124850' LIMIT 3;"
        new QueryHandler()
                .withTarget(result)
                .withStatement(query)
                .perform()

        then:
        result.val[0] == ['acc', 'ERR4796597']
        result.val[0] == ['acc', 'ERR4797168']
        result.val[0] == ['acc', 'ERR4797173']
        result.val == Channel.STOP

    }
}
