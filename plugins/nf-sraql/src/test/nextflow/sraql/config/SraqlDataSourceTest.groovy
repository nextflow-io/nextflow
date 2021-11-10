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
class SraqlDataSourceTest extends Specification {

    def 'should get default config' () {
        given:
        def ds = new SraqlDataSource([:])
        expect:
        ds.source == SraqlDataSource.DEFAULT_SOURCE
    }


    def 'should convert to map' () {
        when:
        def ds = new SraqlDataSource(source:'google-bigquery')
        then:
        ds.toMap().source == 'google-bigquery'
    }

}
