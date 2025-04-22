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

package nextflow.lineage.h2

import nextflow.lineage.model.Annotation
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.DataPath
import nextflow.lineage.model.FileOutput
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.Workflow
import nextflow.lineage.model.WorkflowRun
import nextflow.lineage.config.LineageConfig
import spock.lang.Shared
import spock.lang.Specification

import java.time.OffsetDateTime

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class H2LinStoreTest extends Specification {

    @Shared
    H2LinStore store

    def setupSpec() {
        def uri = "jdbc:h2:mem:testdb;DB_CLOSE_DELAY=-1"
        def config = new LineageConfig([store:[location:uri]])
        store = new H2LinStore().open(config)
    }

    def cleanupSpec() {
        store.close()
    }

    def 'should store and get a value' () {
        given:
        def value = new FileOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "lid://source", "lid://workflow", "lid//task", 1234)
        when:
        store.save('/some/key', value)
        then:
        store.load('/some/key').toString() == value.toString()
    }

    def 'should query' () {
        given:
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def time = OffsetDateTime.now()
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])
        def key2 = "testKey2"
        def value2 = new FileOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "lid://workflow", "lid//task", 1234, time, time, [new Annotation("key1", "value1"), new Annotation("key2", "value2")])
        def key3 = "testKey3"
        def value3 = new FileOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "lid://workflow", "lid//task", 1234, time, time, [new Annotation("key2", "value2"), new Annotation("key3", "value3")])
        def key4 = "testKey4"
        def value4 = new FileOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "lid://workflow", "lid//task", 1234, time, time, [new Annotation("key3", "value3"), new Annotation("key4", "value4")])

        store.save(key, value1)
        store.save(key2, value2)
        store.save(key3, value3)
        store.save(key4, value4)
        when:
        def results = store.search("type=FileOutput&annotations.key=key2&annotations.value=value2")
        then:
        results.size() == 2
    }
}
