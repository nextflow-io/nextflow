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
import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun
import nextflow.data.config.DataConfig
import spock.lang.Shared
import spock.lang.Specification

import java.time.Instant

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
        def value = new DataOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", "cid://run", 1234)
        when:
        store.save('/some/key', value)
        then:
        store.load('/some/key').toString() == value.toString()
    }

    def 'should query' () {
        given:
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow(mainScript, [], "https://nextflow.io/nf-test/", "123456")
        def time = Instant.ofEpochMilli(1234567)
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])
        def key2 = "testKey2"
        def value2 = new DataOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "cid://run", 1234, time, time, [key1: "value1", key2: "value2"])
        def key3 = "testKey3"
        def value3 = new DataOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "cid://run", 1234, time, time, [key2: "value2", key3: "value3"])
        def key4 = "testKey4"
        def value4 = new DataOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "cid://run", 1234, time, time, [key3: "value3", key4: "value4"])

        store.save(key, value1)
        store.save(key2, value2)
        store.save(key3, value3)
        store.save(key4, value4)
        when:
        def results = store.search("type=DataOutput&annotations.key2=value2")
        then:
        results.size() == 2
    }
}
