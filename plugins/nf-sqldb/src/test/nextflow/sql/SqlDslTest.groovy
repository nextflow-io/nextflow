/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.sql

import groovy.sql.Sql
import nextflow.Channel
import nextflow.Session
import nextflow.extension.ChannelExtensionDelegate
import nextflow.plugin.Plugins
import spock.lang.Stepwise
import spock.lang.Timeout
import test.BaseSpec
import test.Dsl2Spec
import test.MockScriptRunner

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class SqlDslTest extends BaseSpec {

    def setup () {
        ChannelExtensionDelegate.reloadExtensionPoints()
    }

    def 'should perform a query and create a channel' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.execute("insert into FOO (id, alpha, omega) values (1, 'hola', 10) ")
        sql.execute("insert into FOO (id, alpha, omega) values (2, 'ciao', 20) ")
        sql.execute("insert into FOO (id, alpha, omega) values (3, 'hello', 30) ")
        and:
        def config = [sql: [db: [test: [url: JDBC_URL]]]]

        when:
        def SCRIPT = '''
            def table = 'FOO'
            def sql = "select * from $table"
            channel.sql.fromQuery(sql, db: "test") 
            '''
        and:
        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()
        then:
        result.val == [1, 'hola', 10]
        result.val == [2, 'ciao', 20]
        result.val == [3, 'hello', 30]
        result.val == Channel.STOP
    }


    def 'should insert channel data into a channel' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(1_000_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        and:
        def config = [sql: [db: [ds1: [url: JDBC_URL]]]]

        when:
        def SCRIPT = '''
            channel
              .of(100,200,300)
              .sqlInsert(into:"FOO", columns:'id', db:"ds1")
            '''
        and:
        def result = new MockScriptRunner(config).setScript(SCRIPT).execute()
        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP
        and:
        def rows =  sql.rows("select id from FOO;")
        and:
        rows.size() == 3
        rows.id == [100, 200, 300]

    }

}
