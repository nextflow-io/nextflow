/*
 * Copyright 2013-2026, Seqera Labs
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
 */

package nextflow.lineage.config


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LineageConfigTest extends Specification {

    def 'should create default config' () {
        when:
        def config = new LineageConfig(Map.of())
        then:
        !config.enabled
        !config.store.location
    }

    def 'should create default with enable' () {
        when:
        def config = new LineageConfig([enabled: true])
        then:
        config.enabled
        !config.store.location
    }

    def 'should create data config with location' () {
        when:
        def config = new LineageConfig(enabled: true, store: [location: "/some/data/store"])
        then:
        config.enabled
        config.store.location == '/some/data/store'
    }
}
