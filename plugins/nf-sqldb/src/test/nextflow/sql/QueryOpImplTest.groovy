package nextflow.sql

import java.nio.file.Files

import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class QueryOpImplTest extends Specification {

    def 'should normalise query' () {
        given:
        def ext = new QueryOpImpl()

        expect:
        ext.normalize('select * from x; ') == 'select * from x;'
        ext.normalize('select * from x ')  == 'select * from x;'
    }

    def 'should create db params' () {
        given:
        def ext = new QueryOpImpl()
        when:
        def props = ext.dbProps([:])
        then:
        props.driver == 'org.h2.Driver'
        props.url == 'jdbc:h2:mem:'

    }

    def 'should connect db' () {
        given:
        def ext = new QueryOpImpl()
        when:
        def conn = ext.connect([:])
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
        new QueryOpImpl()
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
