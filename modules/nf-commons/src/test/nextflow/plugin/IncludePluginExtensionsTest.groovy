package nextflow.plugin

import com.sun.net.httpserver.HttpServer
import groovy.sql.Sql
import nextflow.Channel
import nextflow.extension.ChannelExtensionDelegate
import spock.lang.Specification
import spock.lang.Unroll
import test.Dsl2Spec
import test.MockScriptRunner

import java.nio.file.Files
import java.nio.file.Paths

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class IncludePluginExtensionsTest extends Dsl2Spec {

    def 'should execute custom operator extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()

        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def JDBC_URL = 'jdbc:h2:file:' + folder + '/test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.close()
        def config = [sql: [db: [test: [url: JDBC_URL]]]]

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-sqldb@0.4.0']])

        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()

        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP

        and:
        def rows =  Sql.newInstance(JDBC_URL, 'sa', null).rows("select id from FOO;")
        and:
        rows.size() == 3
        rows.id == [100, 200, 300]

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        ChannelExtensionDelegate.reset()
        Plugins.stop()

        where:
        SCRIPT_TEXT << ['''
            include { sqlInsert } from 'plugin/nf-sqldb'

            channel
              .of(100,200,300)
              .sqlInsert(into:"FOO", columns:'id', db:"test")            
            ''', '''

            include { 
                fromQuery;
                sqlInsert; 
            } from 'plugin/nf-sqldb'

            channel
              .of(100,200,300)
              .sqlInsert(into:"FOO", columns:'id', db:"test")            
            '''
        ]
    }

    def 'should execute custom factory extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def JDBC_URL = 'jdbc:h2:file:' + folder + '/test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.execute("insert into FOO (id, alpha, omega) values (1, 'hola', 10) ")
        sql.execute("insert into FOO (id, alpha, omega) values (2, 'ciao', 20) ")
        sql.execute("insert into FOO (id, alpha, omega) values (3, 'hello', 30) ")
        sql.close()
        def config = [sql: [db: [test: [url: JDBC_URL]]]]
        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-sqldb@0.4.0']])

        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()

        then:
        result
        result.val == ['ID', 'ALPHA', 'OMEGA']
        result.val == [1, 'hola', 10]
        result.val == [2, 'ciao', 20]
        result.val == [3, 'hello', 30]
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionDelegate.reset()

        where:
        SCRIPT_TEXT << ['''
            include { fromQuery } from 'plugin/nf-sqldb'
    
            def table = 'FOO'
            def sql = "select * from $table"            
            channel.fromQuery(sql, db: "test", emitColumns:true)            
            ''','''

            include { fromQuery;  } from 'plugin/nf-sqldb'
            include { sqlInsert } from 'plugin/nf-sqldb'
    
            def table = 'FOO'
            def sql = "select * from $table"            
            channel.fromQuery(sql, db: "test", emitColumns:true)            
            '''
        ]
    }

    def 'should execute custom operator as alias extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def JDBC_URL = 'jdbc:h2:file:' + folder + '/test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.close()

        and:
        def config = [sql: [db: [test: [url: JDBC_URL]]]]

        and:
        def SCRIPT_TEXT = '''
            include { sqlInsert as insertInto } from 'plugin/nf-sqldb'

            channel
              .of(100,200,300)
              .insertInto(into:"FOO", columns:'id', db:"test")            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-sqldb@0.4.0']])

        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()

        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP
        and:
        def rows =  Sql.newInstance(JDBC_URL, 'sa', null).rows("select id from FOO;")
        and:
        rows.size() == 3
        rows.id == [100, 200, 300]

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionDelegate.reset()
    }

    def 'should execute custom factory as alias extension' () {
        given:
        HttpServer server = HttpServer.create(new InetSocketAddress(9900), 0);
        server.createContext("/", new FakeIndexHandler());
        server.start()
        and:
        def folder = Files.createTempDirectory('test')
        Plugins.INSTANCE.mode = 'prod'
        Plugins.INSTANCE.root = folder
        Plugins.INSTANCE.env = [:]
        Plugins.INSTANCE.indexUrl = 'http://localhost:9900/plugins.json'

        and:
        def JDBC_URL = 'jdbc:h2:file:' + folder + '/test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.execute("insert into FOO (id, alpha, omega) values (1, 'hola', 10) ")
        sql.execute("insert into FOO (id, alpha, omega) values (2, 'ciao', 20) ")
        sql.execute("insert into FOO (id, alpha, omega) values (3, 'hello', 30) ")
        sql.close()

        and:
        def config = [sql: [db: [test: [url: JDBC_URL]]]]

        and:
        def SCRIPT_TEXT = '''
            nextflow.enable.dsl=2
            include { fromQuery as sqlQuery } from 'plugin/nf-sqldb'

            def table = 'FOO'
            def sql = "select * from $table"            
            channel.sqlQuery(sql, db: "test", emitColumns:true)            
            '''

        and:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = SCRIPT_TEXT

        when:
        Plugins.setup([plugins: ['nf-sqldb@0.4.0']])

        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()

        then:
        result
        result.val == ['ID', 'ALPHA', 'OMEGA']
        result.val == [1, 'hola', 10]
        result.val == [2, 'ciao', 20]
        result.val == [3, 'hello', 30]
        result.val == Channel.STOP

        cleanup:
        folder?.deleteDir()
        server?.stop(0)
        Plugins.stop()
        ChannelExtensionDelegate.reset()
    }

}
