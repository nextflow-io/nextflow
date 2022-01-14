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
import nextflow.sql.config.SqlDataSource
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class InsertHandlerTest extends Specification {

    def 'should get sql insert from statement' () {
        given:
        def DML = 'insert into foo xyz'
        def handler = new InsertHandler(Mock(SqlDataSource), [statement: DML])
        expect:
        handler.getSqlStatement(null) == DML
    }

    def 'should get sql insert from table and columns' () {
        when:
        def handler1 = new InsertHandler(Mock(SqlDataSource), [into: 'FOO', columns:['x1', 'y2', 'z3']])
        then:
        handler1.getSqlStatement([]) == 'insert into FOO ( x1,y2,z3 ) values ( ?,?,? )'

        when:
        def handler2 = new InsertHandler(Mock(SqlDataSource), [into: 'FOO', columns:'a,b,c'])
        then:
        handler2.getSqlStatement([]) == 'insert into FOO ( a,b,c ) values ( ?,?,? )'
    }

    def 'should fetch columns names from db' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(10_000)
        def TABLE = 'create table FOO(id int primary key, alpha varchar(255), omega int);'
        and:
        def ds = new SqlDataSource([url:JDBC_URL])

        when:
        def handler = new InsertHandler(ds, [into: 'FOO', setup: TABLE])
        then:
        handler.getSqlStatement([]) == 'insert into FOO ( ID,ALPHA,OMEGA ) values ( ?,?,? )'

    }

    def 'should insert data into table' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(10_000)
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        def TABLE = 'create table FOO(id int primary key, alpha varchar(255), omega int);'
        and:
        def ds = new SqlDataSource([url:JDBC_URL])
        and:
        def handler = new InsertHandler(ds, [into: 'FOO', setup: TABLE])

        when:
        handler.perform([1, 'Hello world', 100])
        handler.perform([2, 'Hola mundo', 200])
        handler.perform([3, 'Ciao mondo', 300])
        and:
        handler.close()
        then:
        def result = sql.rows('select * from FOO')
        result.size() == 3
        and:
        result.id == [1,2,3]
        result.alpha == ['Hello world','Hola mundo','Ciao mondo']
        result.omega == [100,200,300]
    }

    def 'should insert single value into table' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(10_000)
        def TABLE = 'create table FOO(id int primary key, alpha varchar(255), omega int);'
        and:
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        def ds = new SqlDataSource([url:JDBC_URL])
        and:
        def handler = new InsertHandler(ds, [into: 'FOO', columns: 'id', setup: TABLE])

        when:
        handler.perform(1)
        handler.perform(2)
        handler.perform(3)
        and:
        handler.close()
        then:
        def result = sql.rows('select id from FOO')
        result.size() == 3
        and:
        result.id == [1,2,3]
    }

    def 'should insert map into table' () {
        given:
        def JDBC_URL = 'jdbc:h2:mem:test_' + Random.newInstance().nextInt(10_000)
        def TABLE = 'create table FOO(id int primary key, alpha varchar(255), omega int);'
        and:
        def sql = Sql.newInstance(JDBC_URL, 'sa', null)
        def ds = new SqlDataSource([url:JDBC_URL])
        and:
        def handler = new InsertHandler(ds, [into: 'FOO', columns: 'id,alpha', setup: TABLE ])

        when:
        handler.perform([id: 1, alpha: 'Hola'])
        handler.perform([id: 2, alpha: 'Ciao'])
        handler.perform([id: 3, alpha: 'Hello'])
        and:
        handler.close()
        then:
        def result = sql.rows('select id, alpha from FOO')
        result.size() == 3
        and:
        result.id == [1,2,3]
        result.alpha == ['Hola', 'Ciao', 'Hello']
    }

}
