package nextflow.sraql


import java.nio.file.Files
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sraql.config.SraqlDataSource
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
        def query = "SELECT *  FROM `nih-sra-datastore.sra.metadata` WHERE organism = 'Mycobacterium tuberculosis' LIMIT 3;"
        new QueryHandler()
                .withTarget(result)
                .withStatement(query)
                .perform()

        then:
        result.val[0] == ['acc', 'ERR3287691']
        result.val[0] == ['acc', 'ERR3287738']
        result.val[0] == ['acc', 'SRR12395057']
        result.val == Channel.STOP

    }
}
