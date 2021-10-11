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

package nextflow.sql.config


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SqlDataSourceTest extends Specification {

    def 'should configure datasource'  () {

        given:
        def CONFIG = '''
            dataSources {
                'default' {
                    url = 'jdbc:h2:mem:' 
                    driver = 'org.h2.Driver'
                    user = 'sa'
                }
                
                myDb1 {
                    url = 'jdbc:h2:mem:' 
                    driver = 'org.h2.Driver'
                    user = 'sa'
                    password = null
                }
            
                myDb2 {
                    url = 'jdbc:postgresql://host:port/database'
                    driver = 'org.postgresql.Driver'
                    user = 'xyz'
                    password = 'foo'                    
                }
                
            }
        '''

    }

    def 'should map url to driver' () {
        given:
        def helper = new SqlDataSource([:])

        expect:
        helper.urlToDriver(JBDC_URL) == DRIVER
        where:
        JBDC_URL                        | DRIVER
        'jdbc:postgresql:database'      | 'org.postgresql.Driver'
        'jdbc:sqlite:database'          | 'org.sqlite.JDBC'
        'jdbc:h2:mem:'                  | 'org.h2.Driver'
        'jdbc:mysql:some-host'          | 'com.mysql.cj.jdbc.Driver'
        'jdbc:mariadb:other-host'       | 'org.mariadb.jdbc.Driver'
        'jdbc:duckdb:'                  | 'org.duckdb.DuckDBDriver'
    }

    def 'should get default config' () {
        given:
        def ds = new SqlDataSource([:])
        expect:
        ds.url == SqlDataSource.DEFAULT_URL
        ds.driver == SqlDataSource.DEFAULT_DRIVER
        ds.user == SqlDataSource.DEFAULT_USER
        ds.password == null
    }


    def 'should get postgresql config' () {
        given:
        def ds = new SqlDataSource([url:'jdbc:postgresql:some-host'])
        expect:
        ds.url == 'jdbc:postgresql:some-host'
        ds.driver == 'org.postgresql.Driver'
        ds.user == SqlDataSource.DEFAULT_USER
        ds.password == null
    }

    def 'should get custom config' () {
        given:
        def config = [
                url:'jdbc:xyz:host-name',
                driver:'this.that.Driver',
                user: 'foo',
                password: 'secret']
        and:
        def ds = new SqlDataSource(config)
        
        expect:
        ds.url == 'jdbc:xyz:host-name'
        ds.driver == 'this.that.Driver'
        ds.user == 'foo'
        ds.password == 'secret'
    }

    def 'should convert to map' () {
        when:
        def ds = new SqlDataSource(url:'x', driver: 'y', user: 'w', password: 'z')
        then:
        ds.toMap().url == 'x'
        ds.toMap().driver == 'y'
        ds.toMap().user == 'w'
        ds.toMap().password == 'z'
    }

    def 'should validate equals & hashcode' () {
        given:
        def ds1 = new SqlDataSource(url:'x', driver: 'y', user: 'w', password: 'z')
        def ds2 = new SqlDataSource(url:'x', driver: 'y', user: 'w', password: 'z')
        def ds3 = new SqlDataSource(url:'p', driver: 'q', user: 'r', password: 'v')
        expect:
        ds1 == ds2
        ds1 != ds3
        and:
        ds1.hashCode() == ds2.hashCode()
        ds1.hashCode() != ds3.hashCode()
    }
}
