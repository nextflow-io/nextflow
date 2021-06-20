package nextflow.sql

import java.nio.file.Files

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sql.config.SqlDatasource
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class QueryHandlerTest extends Specification {

    def 'should normalise query' () {
        given:
        def ext = new QueryHandler()

        expect:
        ext.normalize('select * from x; ') == 'select * from x;'
        ext.normalize('select * from x ')  == 'select * from x;'
    }


    def 'should connect db' () {
        given:
        def ext = new QueryHandler()
        when:
        def conn = ext.connect(new SqlDatasource([:]))
        then:
        conn != null
        cleanup:
        conn.close()
    }

    def 'should perform query' () {
        given:
        def folder = Files.createTempDirectory('test')
        and:
        folder.resolve('test.csv').text  = '''\
        FOO,BAR
        1,hello
        2,ciao
        3,hola
        4,bonjour
        '''.stripIndent()

        when:
        def result = new DataflowQueue()
        def query = "SELECT * FROM CSVREAD('${folder.resolve('test.csv')}') where FOO > 2;"
        new QueryHandler()
                .withTarget(result)
                .withStatement(query)
                .perform()

        then:
        result.val == ['3','hola']
        result.val == ['4','bonjour']
        result.val == Channel.STOP

        cleanup:
        folder.deleteDir()
    }
}
