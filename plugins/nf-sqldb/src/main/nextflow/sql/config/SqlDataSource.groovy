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

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.extension.Bolts

/**
 * Model a dataSource configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
class SqlDataSource {
    public static String DEFAULT_URL = 'jdbc:h2:mem:'
    public static String DEFAULT_DRIVER = 'org.h2.Driver'
    public static String DEFAULT_USER = 'sa'

    static SqlDataSource DEFAULT = new SqlDataSource(Collections.emptyMap())

    String driver
    String url
    String user
    String password

    SqlDataSource(Map opts) {
        this.url = opts.url ?: DEFAULT_URL
        this.driver = opts.driver ?: urlToDriver(url) ?: DEFAULT_DRIVER
        this.user = opts.user ?: DEFAULT_USER
        this.password = opts.password
    }

    SqlDataSource(Map opts, SqlDataSource fallback) {
        this.url = opts.url ?: fallback.url ?: DEFAULT_URL
        this.driver = opts.driver ?: urlToDriver(url) ?: fallback.driver ?: DEFAULT_DRIVER
        this.user = opts.user ?: fallback.user ?: DEFAULT_USER
        this.password = opts.password ?: fallback.password
    }


    protected String urlToDriver(String url) {
        if( !url ) return null
        if( !url.startsWith('jdbc:') ) throw new IllegalArgumentException("Invalid database JDBC connection url: $url")
        switch (url.tokenize(':')[1]) {
            case 'h2': return 'org.h2.Driver'
            case 'sqlite': return 'org.sqlite.JDBC'
            case 'mysql': return 'com.mysql.cj.jdbc.Driver'
            case 'mariadb': return 'org.mariadb.jdbc.Driver'
            case 'postgresql': return 'org.postgresql.Driver'
            case 'duckdb': return 'org.duckdb.DuckDBDriver'
        }
        return null
    }

    Map toMap() {
        final result = new HashMap(10)
        if( url )
            result.url = url
        if( driver )
            result.driver = driver
        if( user || password ) {
            result.user = user
            result.password = password
        }
        return result
    }

    @Override
    String toString() {
        return "SqlDataSource[url=$url; driver=$driver; user=$user; password=${Bolts.redact(password)}]"
    }
}
