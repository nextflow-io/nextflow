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

package nextflow.sraql.config

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SraqlConfigTest extends Specification {

    def 'should create empty config' () {
        when:
        def config = new SraqlConfig(null)
        then:
        config.getDataSource('default') == SraqlDataSource.DEFAULT
    }

    def 'should create custom config' () {
        given:
        def map = [myDatabase: [
                    url: 'X1',
                    driver: 'X2',
                    user: 'X3',
                    password: 'X4'
                ]]
        when:
        def config = new SraqlConfig(map)
        then:
        with(config.getDataSource('myDatabase')) {
            url == 'X1'
            driver == 'X2'
            user == 'X3'
            password == 'X4'
        }
    }

    def 'should create config from default' () {
        given:
        def map = [myDatabase: [
                url: 'jdbc:postgresql:host.name',
        ]]
        when:
        def config = new SraqlConfig(map)
        then:
        with(config.getDataSource('myDatabase')) {
            url == 'jdbc:postgresql:host.name'
            driver == 'org.postgresql.Driver'
            user == SqlDataSource.DEFAULT_USER
            password == null
        }
        and:
        with(config.getDataSource('default')) {
            url == SqlDataSource.DEFAULT_URL
            driver == SqlDataSource.DEFAULT_DRIVER
            user == SqlDataSource.DEFAULT_USER
            password == null
        }
    }

    def 'should override default config' () {
        given:
        def map = [
                'default': [url:'jdbc:foo:mem', driver: 'org.foo.Driver', user: 'user-x', password: 'pass-y'],
                myDatabase: [:]]
        when:
        def config = new SraqlConfig(map)
        then:
        with(config.getDataSource('myDatabase')) {
            url == 'jdbc:foo:mem'
            driver == 'org.foo.Driver'
            user == 'user-x'
            password == 'pass-y'
        }
        and:
        with(config.getDataSource('default')) {
            url == 'jdbc:foo:mem'
            driver == 'org.foo.Driver'
            user == 'user-x'
            password == 'pass-y'
        }
    }

    def 'should override default config/2' () {
        given:
        def map = [
                'default': [url:'jdbc:foo:mem', driver: 'org.foo.Driver', user: 'user-x', password: 'pass-y'],
                myDatabase: [url:'jdbc:foo:mem:custom']]
        when:
        def config = new SraqlConfig(map)
        then:
        with(config.getDataSource('myDatabase')) {
            url == 'jdbc:foo:mem:custom'
            driver == 'org.foo.Driver'
            user == 'user-x'
            password == 'pass-y'
        }
        and:
        and:
        with(config.getDataSource('default')) {
            url == 'jdbc:foo:mem'
            driver == 'org.foo.Driver'
            user == 'user-x'
            password == 'pass-y'
        }
    }


    def 'should validate equals & hashCode' () {
        given:
        def foo = [ 'myDatabase': [url:'jdbc:foo:mem', driver: 'org.foo.Driver', user: 'user-x', password: 'pass-y'] ]
        def bar = [ 'myDatabase': [url:'jdbc:bar:mem', driver: 'org.bar.Driver', user: 'user-x', password: 'pass-y'] ]
        and:
        def cfg1 = new SraqlConfig(foo)
        def cfg2 = new SraqlConfig(foo)
        def cfg3 = new SraqlConfig(bar)

        expect:
        cfg1 == cfg2
        cfg1 != cfg3
        and:
        cfg1.hashCode() == cfg2.hashCode()
        cfg1.hashCode() != cfg3.hashCode()


    }
}
