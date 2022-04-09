package nextflow.sql

import java.nio.file.Files

import groovy.sql.Sql
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.sql.config.SqlDataSource
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
        def conn = ext.connect(new SqlDataSource([:]))
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

    def 'should emit header when perform query' () {
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
                .withOpts(emitColumns: true)
                .perform()

        then:
        result.val == ['FOO','BAR']
        result.val == ['3','hola']
        result.val == ['4','bonjour']
        result.val == Channel.STOP

        cleanup:
        folder.deleteDir()
    }
    
    def 'should append limit and offset' () {
        given:
        def ext = new QueryHandler()
        
        expect:
        ext.makePaginationStm('select * from FOO')  == 'select * from FOO LIMIT ? OFFSET ?;'
        ext.makePaginationStm('select * from FOO  ;  ')  == 'select * from FOO LIMIT ? OFFSET ?;'

        when:
        ext.makePaginationStm('select * from offset')
        then:
        thrown(IllegalArgumentException)

        when:
        ext.makePaginationStm('select * from limit')
        then:
        thrown(IllegalArgumentException)
    }

    def 'should test paginated query' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(10_000)
        def TABLE = 'create table FOO(id int primary key, alpha varchar(255));'
        def ds = new SqlDataSource([url:JDBC_URL])
        and:
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        sql.execute(TABLE)
        for( int x : 1..13 ) {
            def params = [x, "Hello $x".toString()]
            sql.execute("insert into FOO (id, alpha) values (?,?);", params)
        }

        when:
        def result = new DataflowQueue()
        def query = "SELECT id, alpha from FOO order by id "
        def handler = new QueryHandler()
                .withTarget(result)
                .withStatement(query)
                .withDataSource(ds)
                .withOpts(batchSize: 5)
                .perform()
        then:
        handler.batchSize() == 5
        handler.queryCount() == 3
        and:
        result.length() == 14  // <-- 13 + the stop signal value
        and:
        result.getVal() == [1, 'Hello 1']
        result.getVal() == [2, 'Hello 2']
        result.getVal() == [3, 'Hello 3']
        result.getVal() == [4, 'Hello 4']
        result.getVal() == [5, 'Hello 5']
        result.getVal() == [6, 'Hello 6']
        result.getVal() == [7, 'Hello 7']
        result.getVal() == [8, 'Hello 8']
        result.getVal() == [9, 'Hello 9']
        result.getVal() == [10, 'Hello 10']
        result.getVal() == [11, 'Hello 11']
        result.getVal() == [12, 'Hello 12']
        result.getVal() == [13, 'Hello 13']
        result.getVal() == Channel.STOP

        when:
        def result2 = new DataflowQueue()
        def query2 = "SELECT id, alpha from FOO order by id "
        new QueryHandler()
                .withTarget(result2)
                .withStatement(query2)
                .withDataSource(ds)
                .withOpts(batchSize: 5, emitColumns: true)
                .perform()
        then:
        result2.length() == 15  // <-- 13 + columns name tuple + the stop signal value
        and:
        result2.getVal() == ['ID', 'ALPHA']
        result2.getVal() == [1, 'Hello 1']
        result2.getVal() == [2, 'Hello 2']
    }
}
