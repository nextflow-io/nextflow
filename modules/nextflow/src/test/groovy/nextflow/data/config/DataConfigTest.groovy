/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.data.config

import java.nio.file.Path

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataConfigTest extends Specification {

    def 'should create default config' () {
        when:
        def config = new DataConfig(Map.of())
        then:
        config.store.location == Path.of('.').resolve('data').toAbsolutePath().normalize()
        config.store.logLocation == null
        !config.enabled
    }

    def 'should create default with enable' () {
        when:
        def config = new DataConfig([enabled: true])
        then:
        config.store.location == Path.of('.').resolve('data').toAbsolutePath().normalize()
        config.store.logLocation == null
        config.enabled
    }

    def 'should create data config with location' () {
        when:
        def config = new DataConfig(enabled: true, store: [location: "/some/data/store", logLocation: "/some/data/.history"])
        then:
        config.store.location == Path.of("/some/data/store")
        config.store.logLocation == Path.of("/some/data/.history")
        config.enabled
    }
}
