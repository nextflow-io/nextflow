/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.data.cid.h2

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.Output
import nextflow.data.cid.model.DataType
import nextflow.data.config.DataConfig
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class H2CidStoreTest extends Specification {

    @Shared
    H2CidStore store

    def setupSpec() {
        def uri = "jdbc:h2:mem:testdb;DB_CLOSE_DELAY=-1"
        def config = new DataConfig([store:[location:uri]])
        store = new H2CidStore().open(config)
    }

    def cleanupSpec() {
        store.close()
    }

    def 'should store and get a value' () {
        given:
        def value = new Output(DataType.TaskOutput, "/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", 1234)
        when:
        store.save('/some/key', value)
        then:
        store.load('/some/key').toString() == value.toString()
    }

}
