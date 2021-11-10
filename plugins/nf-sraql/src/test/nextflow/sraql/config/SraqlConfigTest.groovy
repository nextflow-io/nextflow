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
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
class SraqlConfigTest extends Specification {

    def 'should create empty config' () {
        when:
        def config = new SraqlConfig(null)
        then:
        println(config)
    }

    def 'should create config from default' () {
        given:
        def map = [google: [
                source: 'google-bigquery',
        ]]

        when:
        def config = new SraqlConfig(map)

        then:
        println("---\n${config}\n---")
//        with(config.getDataSource('google')) {
//            source == 'google-bigquery'
//        }
//
//        and:
//        with(config.getDataSource('default')) {
//            source == 'google-bigquery'
//        }
    }

    def 'should override default config' () {
        given:
        def map = [
                'default': [source:'google-bigquery'],
                google: [:]]
        when:
        def config = new SraqlConfig(map)

        then:
        with(config.getDataSource('myDatabase')) {
            source == 'google-bigquery'
        }

        and:
        with(config.getDataSource('default')) {
            source == 'google-bigquery'
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
