/*
 * Copyright 2020-2021, Seqera Labs
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

package nextflow.sraql

import groovy.sql.Sql
import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel
import nextflow.Session
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(10)
class ChannelSraqlExtensionTest extends Specification {

    def 'should create channel from query' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(1_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        and:
        sql.execute('create table FOO(id int primary key, alpha varchar(255), omega int);')
        sql.execute("insert into FOO (id, alpha, omega) values (1, 'hola', 10) ")
        sql.execute("insert into FOO (id, alpha, omega) values (2, 'ciao', 20) ")
        sql.execute("insert into FOO (id, alpha, omega) values (3, 'hello', 30) ")
        and:
        def session = Mock(Session) {
            getConfig() >> [sql: [db: [test: [url: JDBC_URL]]]]
        }
        def sqlExtension = new ChannelSraqlExtension(); sqlExtension.init(session)

        when:
        def result = sqlExtension.fromQuery('select * from FOO', db: 'test')
        then:
        result.val == [1, 'hola', 10]
        result.val == [2, 'ciao', 20]
        result.val == [3, 'hello', 30]
        result.val == Channel.STOP

        when:
        result = sqlExtension.fromQuery('select alpha, omega from FOO where id=3', db: 'test')
        then:
        result.val == ['hello', 30]
        result.val == Channel.STOP
    }
}
